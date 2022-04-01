using ManagedBass;
using ManagedBass.Fx;
using osu.Framework.Audio.Track;
using System.Diagnostics;
using System.Runtime.InteropServices;

Bass.Init();

/*
 * 
 * references links:
 * 
 * ITU-R BS.1770-4, Algorithms to measure audio programme loudness and true-peak audio level (https://www.itu.int/rec/R-REC-BS.1770-4-201510-I)
 * EBU R 128, Loudness normalisation and permitted maximum level of audio signals (https://tech.ebu.ch/publications/e128)
 * EBU TECH 3341, Loudness Metering: ‘EBU Mode’ metering to supplement loudness normalisation in accordance with EBU R 128 (https://tech.ebu.ch/publications/tech3341)
 * 
 * 
 */

//needed to apply effects
var version = BassFx.Version;

//absolute silence in ITU_R_BS._1770._4 recommendation (LUFS)
const double ABSOLUTE_SILENCE = -70;

//same reference level as apple music (LUFS)
const double REFERENCE_LEVEL = -16;

//sampling frequency at which the filter coeffs are provided from the ITU_R_BS._1770._4 recommendation
const double REFERENCE_SAMPLING_FREQUENCY = 48000.0;

Console.Write("Write the path to the audio: ");
string path = Console.ReadLine();
Console.WriteLine("\n");

if (path == "")
{
    Console.WriteLine("Invalid path");
    Console.ReadKey();
    Environment.Exit(0);

}
else
{
    path = path.Replace('\\', '/');
}

Stopwatch stopWatch = new();
stopWatch.Start();
/*
 * 
 * Load track into Bass and initialize majority of variables needed
 * 
 */
//load the track and read it's info
var decodeStream = Bass.CreateStream(path, 0, 0, BassFlags.Decode | BassFlags.Float);

if (decodeStream == 0)
{
    Console.WriteLine("Error creating the audio stream: " + Bass.LastError);
    Console.ReadKey();
    Environment.Exit(0);
}

Bass.ChannelGetInfo(decodeStream, out ChannelInfo info);
long length;

Console.WriteLine("Sample frequency: " + info.Frequency + "\n");

//100 ms window
int samplesPerWindow = (int)(info.Frequency * 0.1f * info.Channels);
int bytesPerWindow = samplesPerWindow * TrackBass.BYTES_PER_SAMPLE;

//create a 100ms buffer and read the first segment of the track
float[] sampleBuffer = new float[samplesPerWindow];
length = Bass.ChannelGetData(decodeStream, sampleBuffer, bytesPerWindow);

//number of samples in a 400ms window for 1 channel
int totalWindowLength = samplesPerWindow * 4 / info.Channels;

List<double>[] squaredSegments = new List<double>[info.Channels].Select(item => new List<double> { }).ToArray();


double peakAmp = 0;
//for thread safety when writing peakAmp
object maxLock = new object();

/*
 * 
 *  this section allocates the coefficients fot the “K” frequency weighting
 *  requested from the from ITU_R_BS._1770._4 recommendation.
 * 
 */

//pre-filter coeffs to model a spherical head
double[] headA;
double[] headB;

//high-pass coeffs
double[] highPassA;
double[] highPassB;

if (info.Frequency == 48000)
{
    //filter coeffs for 48000hz sampling rate (from ITU_R_BS._1770._4 guidelines)
    headA = new double[] { -1.69065929318241, 0.73248077421585 };
    headB = new double[] { 1.53512485958697, -2.69169618940638, 1.19839281085285 };

    highPassA = new double[] { -1.99004745483398, 0.99007225036621 };
    highPassB = new double[] { 1.0, -2.0, 1.0 };
}
else if (info.Frequency == 44100)
{
    //filter coeffs for 44100hz sampling rate precalculated for speed
    headA = new double[] { -1.6636551132560202, 0.7125954280732254};
    headB = new double[] { 1.5308412300503478, -2.650979995154729, 1.1690790799215869 };

    highPassA = new double[] { -1.9891696736297957, 0.9891990357870394 };
    highPassB = new double[] { 0.9995600645425144, -1.999120129085029, 0.9995600645425144 };
}
else
{

    //for the rest of the sampling rates coeffs are recalculated at runtime
    headA = new double[] { -1.69065929318241, 0.73248077421585 };
    headB = new double[] { 1.53512485958697, -2.69169618940638, 1.19839281085285 };

    highPassA = new double[] { -1.99004745483398, 0.99007225036621 };
    highPassB = new double[] { 1.0, -2.0, 1.0 };

    filterCoeffsRecalc(info.Frequency, headA, headB);
    filterCoeffsRecalc(info.Frequency, highPassA, highPassB);
}


/*
 * 
 * Multi-threaded calculus of the sum of squared samples for each channel after the
 * “K” frequency weighting is applied.
 * Peak amplitude is also found in a naive way by storing the max absolute value
 * of the samples.
 * 
 */

//list of started tasks
List<Task> squareSegmentTasks = new List<Task>();
int segNumber = 0;

//read the full track excluding last segment if it does not fill the buffer
while (length == bytesPerWindow)
{
    //start a task for every 100ms window
    for (int i = 0; i < info.Channels; i++)
    {
        //task parameters
        int currChannel = i;
        int currSegment = segNumber;
        float[] currentBuffer = new float[sampleBuffer.Length];
        Array.Copy(sampleBuffer, currentBuffer, sampleBuffer.Length);

        //pre allocation of the memory that the task will write into
        squaredSegments[i].Add(0);
        squareSegmentTasks.Add( Task.Run(() => segmentSquaredByChannel(currChannel, currentBuffer, currSegment)));

    }

    //read next segment
    length = Bass.ChannelGetData(decodeStream, sampleBuffer, bytesPerWindow);
    segNumber++;

}
Bass.StreamFree(decodeStream);


/*
 * 
 * calculation of the mean square for each channel over 400ms windows overlapping by 75%
 * as per ITU_R_BS._1770._4 recommendation.
 * 
 */

//list of the squared mean for the overlapping 400ms windows
List<double>[] squaredMeanByChannel = new List<double>[info.Channels].Select(item => new List<double> { }).ToArray();

try
{
    Task.WaitAll(squareSegmentTasks.ToArray());
}
catch (AggregateException ae)
{
    Console.WriteLine("One or more exceptions occurred: ");
    foreach (var ex in ae.Flatten().InnerExceptions)
        Console.WriteLine("   {0}", ex.Message);
}

if (squaredSegments[0].Count == 0)
{
    Console.WriteLine("The file is empty");
    Environment.Exit(0);
}

List<Task> squaremeanTasks = new List<Task>();

//start from 400ms in since windowedSquaredMean reads previous segments
for (int i = 3; i < squaredSegments[0].Count; i++)
{
    for(int j = 0; j < info.Channels; j++)
    {
        //task paramenters
        int currChannel = j;
        int currSegment = i;

        //pre allocation of the memory that the task will write into
        squaredMeanByChannel[currChannel].Add(0);
        squaremeanTasks.Add(Task.Run(() => windowedSquaredMean(currChannel, currSegment)));
    }
}


/*
 * 
 * Pre gating loudness of each 400ms window in LUFS as per ITU_R_BS._1770._4 recommendation.
 * this specific implementation will not work well for channel counts above 5.
 * 
 */

//loudness of each 400ms window when all channels are summed
List<double> blockLoudness = new List<double> { };

try
{
    Task.WaitAll(squaremeanTasks.ToArray());
}
catch (AggregateException ae)
{
    Console.WriteLine("One or more exceptions occurred: ");
    foreach (var ex in ae.Flatten().InnerExceptions)
        Console.WriteLine("   {0}", ex.Message);
}

if (info.Channels > 3 && info.Channels < 6)
{
    for (int i = 0; i < squaredMeanByChannel[0].Count; i++)
    {
        double tempSum = 0;

        for (int j = 0; j < info.Channels - 2; j++)
            tempSum += squaredMeanByChannel[j][i];
        for (int j = info.Channels - 2; j < info.Channels - 2; j++)
            tempSum += squaredMeanByChannel[j][i] * 1.41;

        blockLoudness.Add(-0.691 + 10 * Math.Log10(tempSum));
    }
}
else if(info.Channels >= 6)
{
    for (int i = 0; i < squaredMeanByChannel[0].Count; i++)
    {
        double tempSum = 0;

        for (int j = 0; j < info.Channels; j++)
        {
            if (j == 4 || j == 5)
                tempSum += squaredMeanByChannel[j][i] * 1.41;
            else
                tempSum += squaredMeanByChannel[j][i];
        }

        blockLoudness.Add(-0.691 + 10 * Math.Log10(tempSum));
    }
}
else
{
    for (int i = 0; i < squaredMeanByChannel[0].Count; i++)
    {
        double tempSum = 0;

        for (int j = 0; j < info.Channels; j++)
            tempSum += squaredMeanByChannel[j][i];

        blockLoudness.Add(-0.691 + 10 * Math.Log10(tempSum));
    }
}

/*
 * 
 * gated loudness calc as per ITU_R_BS._1770._4 recommendation.
 * 
 */

double relativeGate = relativeGateCalc();

if (squaredMeanByChannel[0].Count == 0)
{
    Console.WriteLine("The file is completely silent");
    Console.ReadKey();
    Environment.Exit(0);
}

double currLoudness = gatedLoudnessCalc(relativeGate);

stopWatch.Stop();

TimeSpan ts = stopWatch.Elapsed;
string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
    ts.Hours, ts.Minutes, ts.Seconds,
    ts.Milliseconds);
Console.WriteLine("RunTime(H:M:S:Ms) " + elapsedTime + "\n");


Console.WriteLine("Original loundess(LUFS): " + currLoudness);

//Gain in LUFS that needs to be applied to the track
double gain = REFERENCE_LEVEL - currLoudness;

Console.WriteLine("Gain needed: " + gain);
Console.WriteLine("Peak amplitude: "+ peakAmp);







//create a new stream for testing
var decodeStream2 = Bass.CreateStream(path);


//add a Volume effect to the stream (Managedbass has the wrong implementation of this so a hard cast to 9 is needed)
int fx = Bass.ChannelSetFX(decodeStream2, (EffectType)9, 0);

//create an instance of the param structure for volume and create a pointer to it
BASS_FX_VOLUME_PARAM newParams = new BASS_FX_VOLUME_PARAM();
IntPtr re = Marshal.AllocHGlobal(Marshal.SizeOf(newParams));

//set the wanted Volume parameters
newParams.fCurrent = 1;
newParams.fTarget = (float)Math.Pow(10, gain / 20); //inverse of loudness calculation
newParams.fTime = 0;
Marshal.StructureToPtr(newParams, re, true);

//Update the volume effect with the new params
Bass.FXSetParameters(fx, re);

//if the peak amplitude after the gain would clip apply a compressor
if (Math.Pow(10, gain  / 20) * peakAmp >= 1)
{

    Console.WriteLine("Applying a compressor to prevent clipping");
    int compressor = Bass.ChannelSetFX(decodeStream2, EffectType.Compressor, 1);
    CompressorParameters compParams = new();
    Bass.FXGetParameters(compressor, compParams);
    compParams.fAttack = 0.01f;
    compParams.fGain = 0;
    compParams.fRatio = 30f;
    compParams.fRelease = 200;
    compParams.fThreshold = -18;
    Bass.FXSetParameters(compressor, compParams);

}

Console.WriteLine("Volume multiplier: " + newParams.fTarget + "\n\n");

Bass.ChannelPlay(decodeStream2);
Console.WriteLine("Press a key to stop the audio");

Console.ReadKey();
Bass.StreamFree(decodeStream2);


//for a window apply pre filters and find it's sum of all samples + update peak amplitude value if needed
void segmentSquaredByChannel(int channel, float[] data, int segmentIndex)
{
    //Variables to apply the 1st pre-filter
    double pastX0 = 0;
    double pastX1 = 0;

    double pastZ0 = 0;
    double pastZ1 = 0;


    //Variables for the high-pass filter
    double pastZlow0 = 0;
    double pastZlow1 = 0;

    double pastY0 = 0;
    double pastY1 = 0;


    double partialSample = 0;
    double localMax = 0;

    for (int s = channel; s < data.Length; s+= info.Channels)
    {
        /*
         * 
         * “K” frequency weighting for each sample
         * 
         */

        //apply the 1st pre-filter to the sample
        double yuleSample = headB[0] * data[s] + headB[1] * pastX0 + headB[2] * pastX1
            - headA[0] * pastZ0 - headA[1] * pastZ1;

        pastX1 = pastX0;
        pastZ1 = pastZ0;
        pastX0 = data[s];
        pastZ0 = yuleSample;

        //apply the high-pass filter to the sample
        double tempsample = highPassB[0] * yuleSample + highPassB[1] * pastZlow0 + highPassB[2] * pastZlow1
            - highPassA[0] * pastY0 - highPassA[1] * pastY1;

        pastZlow1 = pastZlow0;
        pastY1 = pastY0;
        pastZlow0 = yuleSample;
        pastY0 = tempsample;

        /*
         * sum of squared samples and localMax update
         */
        partialSample += tempsample * tempsample;

        if (Math.Abs(data[s]) > localMax)
        {
            localMax = Math.Abs(data[s]);
        }

    }

    //memorize the sum of squared samples for the given data
    squaredSegments[channel][segmentIndex] = partialSample;

    //thread safe update of peak Amplitude
    if (localMax > peakAmp)
    {
        lock (maxLock)
        {
            if (localMax > peakAmp)
            {
                peakAmp = localMax;
            }

        }

    }

}

//squared mean of a 400ms segment
void windowedSquaredMean(int channel, int segmentIndex)
{
    squaredMeanByChannel[channel][segmentIndex - 3] = (squaredSegments[channel][segmentIndex - 3] 
        + squaredSegments[channel][segmentIndex - 2] 
        + squaredSegments[channel][segmentIndex - 1] 
        + squaredSegments[channel][segmentIndex]) / totalWindowLength;
}

//calc of the relative gate for the loudness as per ITU_R_BS._1770._4 recommendation
double relativeGateCalc()
{

    double tempTotLoudness = 0;
    int nonSilenceSegments = 0;
    /*
     * removal of segments below the silence threshold
     */
    for (int i = 0; i < blockLoudness.Count; i++)
    {
        if (blockLoudness[i] > ABSOLUTE_SILENCE)
        {
            for (int j = 0; j < info.Channels; j++)
            {
                tempTotLoudness += squaredMeanByChannel[j][i];
            }
            nonSilenceSegments++;
        }
        else
        {
            blockLoudness.RemoveAt(i);
            for (int j = 0; j < info.Channels; j++)
            {
                squaredMeanByChannel[j].RemoveAt(i);
            }
        }
    }

    //relative gate as per ITU_R_BS._1770._4 recommendation
    return -0.691 + 10 * Math.Log10(tempTotLoudness / nonSilenceSegments) - 10;
}

//calc of the double gated loudness as per ITU_R_BS._1770._4 recommendation
double gatedLoudnessCalc(double relativeGate)
{

    double tempTotLoudness = 0;
    int aboveGatesSegments = 0;

    /*
     * removal of segments below the relative gate threshold
     */
    for (int i = 0; i < blockLoudness.Count; i++)
    {
        if (blockLoudness[i] > relativeGate)
        {
            for (int j = 0; j < info.Channels; j++)
            {
                tempTotLoudness += squaredMeanByChannel[j][i];
            }
            aboveGatesSegments++;
        }
    }

    return -0.691 + 10 * Math.Log10(tempTotLoudness / aboveGatesSegments);
}

//calc of filter coeffs for sampling rates different from 48000
void filterCoeffsRecalc(double samplingFreq, double[] a, double[] b)
{

    /*
     * adapted from https://github.com/klangfreund/LUFSMeter/blob/master/filters/SecondOrderIIRFilter.cpp
     * derivation of the equations can be found at https://github.com/klangfreund/LUFSMeter/blob/master/docs/developmentNotes/111222_filter_coefficients/111222_my_notes_to_the_calculation_of_the_filter_coefficients.tif
     *
     */
    double KoverQ = (2.0 - 2.0 * a[1]) / (a[1] - a[0] + 1.0);
    double K = Math.Sqrt((a[0] + a[1] + 1.0) / (a[1] - a[0] + 1.0));
    double Q = K / KoverQ;
    double ArctanK = Math.Atan(K);
    double VB = (b[0] - b[2]) / (1.0 - a[1]);
    double VH = (b[0] - b[1] + b[2]) / (a[1] - a[0] + 1.0);
    double VL = (b[0] + b[1] + b[2]) / (a[0] + a[1] + 1.0);

    double newK = Math.Tan(ArctanK * REFERENCE_SAMPLING_FREQUENCY / samplingFreq);
    double commonFactor = 1.0 / (1.0 + newK / Q + newK * newK);
    
    b[0] = (VH + VB * newK / Q + VL * newK * newK) * commonFactor;
    b[1] = 2.0 * (VL * newK * newK - VH) * commonFactor;
    b[2] = (VH - VB * newK / Q + VL * newK * newK) * commonFactor;
    a[0] = 2.0 * (newK * newK - 1.0) * commonFactor;
    a[1] = (1.0 - newK / Q + newK * newK) * commonFactor;

}


//struct for the Volume FX parameters as per Bass specification
struct BASS_FX_VOLUME_PARAM
{
    public BASS_FX_VOLUME_PARAM() { }
    public float fTarget = 1; //new volume to reach
    public float fCurrent = 1;  //current volume
    public float fTime = 0; //time to reach fTarget
    public int lCurve = 0;  //curve used to reach fTarget
};
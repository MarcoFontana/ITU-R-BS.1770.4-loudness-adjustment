using ManagedBass;
using ManagedBass.Fx;
using osu.Framework.Audio.Track;
using System.Diagnostics;
using System.Runtime.InteropServices;

Bass.Init();

//needed to apply effects
var version = BassFx.Version;

//absolute silence in ITU_R_BS._1770._4 guidelines (LUFS)
const double ABSOLUTE_SILENCE = -70;

//same reference level as apple music (LUFS)
const double REFERENCE_LEVEL = -16;

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


//pre-filter coeffs to model a spherical head
double[] headA;
double[] headB;


//high-pass coeffs
double[] highPassA;
double[] highPassB;

headA = new double[] { -1.69065929318241, 0.73248077421585 };
headB = new double[]{ 1.53512485958697, -2.69169618940638, 1.19839281085285 };

highPassA = new double[] { -1.99004745483398, 0.99007225036621 };
highPassB = new double[] { 1.0, -2.0, 1.0 };

//list of started tasks
List<Task> squareSegmentTasks = new List<Task>();
int segNumber = 0;

//read the full track excluding last segment if it does not fill the buffer
while (length == bytesPerWindow)
{
    //start a task for every 100ms window
    for (int i = 0; i < info.Channels; i++)
    {

        int currChannel = i;
        int currSegment = segNumber;
        float[] currentBuffer = new float[sampleBuffer.Length];
        Array.Copy(sampleBuffer, currentBuffer, sampleBuffer.Length);
        squaredSegments[i].Add(0);
        squareSegmentTasks.Add( Task.Run(() => segmentSquaredByChannel(currChannel, currentBuffer, currSegment)));

    }

    //read next segment
    length = Bass.ChannelGetData(decodeStream, sampleBuffer, bytesPerWindow);
    segNumber++;

}
Bass.StreamFree(decodeStream);

//list of the squared mean for the overlapping 400ms windows
List<double>[] squaredMeanByChannel = new List<double>[info.Channels].Select(item => new List<double> { }).ToArray();

Task.WaitAll(squareSegmentTasks.ToArray());

if (squaredSegments[0].Count == 0)
{
    Console.WriteLine("The file is empty");
    Environment.Exit(0);
}

List<Task> squaremeanTasks = new List<Task>();

for (int i = 3; i < squaredSegments[0].Count; i++)
{
    for(int j = 0; j < info.Channels; j++)
    {
        int currChannel = j;
        int currSegment = i;
        squaredMeanByChannel[currChannel].Add(0);
        squaremeanTasks.Add(Task.Run(() => windowedSquaredMean(currChannel, currSegment)));
    }
}

//loudness of each 400ms window when all channels are summed
List<double> blockLoudness = new List<double> { };

Task.WaitAll(squaremeanTasks.ToArray());

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
    compParams.fRatio = 1.5f;
    compParams.fRelease = 200;
    compParams.fThreshold = -50;
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


        partialSample += tempsample * tempsample;

        if (Math.Abs(data[s]) > localMax)
        {
            localMax = Math.Abs(data[s]);
        }

    }

    squaredSegments[channel][segmentIndex] = partialSample;

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

//calc of the relative gate for the loudness as per ITU_R_BS._1770._4 guidelines
double relativeGateCalc()
{

    double tempTotLoudness = 0;
    int nonSilenceSegments = 0;

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

    return -0.691 + 10 * Math.Log10(tempTotLoudness / nonSilenceSegments) - 10;
}

//calc of the double gated loudness as per ITU_R_BS._1770._4 guidelines
double gatedLoudnessCalc(double relativeGate)
{

    double tempTotLoudness = 0;
    int aboveGatesSegments = 0;

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


//struct for the Volume FX parameters as per Bass specification
struct BASS_FX_VOLUME_PARAM
{
    public BASS_FX_VOLUME_PARAM() { }
    public float fTarget = 1; //new volume to reach
    public float fCurrent = 1;  //current volume
    public float fTime = 0; //time to reach fTarget
    public int lCurve = 0;  //curve used to reach fTarget
};
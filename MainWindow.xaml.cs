using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows;
using NAudio.Wave;
using OxyPlot;
using OxyPlot.Series;
using MathNet.Numerics;
using MathNet.Numerics.IntegralTransforms;

namespace FFTDemoApp
{
    /// <summary>
    /// MainWindow.xaml 的交互逻辑
    /// </summary>
    public partial class MainWindow : System.Windows.Window
    {

        #region 变量定义

        int size = 4096;
        int fs = 4000 * 10;
        double[] reFFT;
        double[] imFFT;

        // FFT
        public static void FFT(double[] inputRe, double[] inputIm, out double[] outputRe, out double[] outputIm, int bitSize)
        {
            int dataSize = 1 << bitSize;
            int[] reverseBitArray = BitScrollArray(dataSize);

            outputRe = new double[dataSize];
            outputIm = new double[dataSize];

            for (int i = 0; i < dataSize; i++)
            {
                outputRe[i] = inputRe[reverseBitArray[i]];
                outputIm[i] = inputIm[reverseBitArray[i]];
            }

            for (int stage = 1; stage <= bitSize; stage++)
            {
                int butterflyDistance = 1 << stage;
                int numType = butterflyDistance >> 1;
                int butterflySize = butterflyDistance >> 1;

                double wRe = 1.0;
                double wIm = 0.0;
                double uRe = Math.Cos(System.Math.PI / butterflySize);
                double uIm = -Math.Sin(System.Math.PI / butterflySize);

                for (int type = 0; type < numType; type++)
                {
                    for (int j = type; j < dataSize; j += butterflyDistance)
                    {
                        int jp = j + butterflySize;
                        double tempRe = outputRe[jp] * wRe - outputIm[jp] * wIm;
                        double tempIm = outputRe[jp] * wIm + outputIm[jp] * wRe;
                        outputRe[jp] = outputRe[j] - tempRe;
                        outputIm[jp] = outputIm[j] - tempIm;
                        outputRe[j] += tempRe;
                        outputIm[j] += tempIm;
                    }
                    double tempWRe = wRe * uRe - wIm * uIm;
                    double tempWIm = wRe * uIm + wIm * uRe;
                    wRe = tempWRe;
                    wIm = tempWIm;
                }
            }
        }

        private static int[] BitScrollArray(int arraySize)
        {
            int[] reBitArray = new int[arraySize];
            int arraySizeHarf = arraySize >> 1;

            reBitArray[0] = 0;
            for (int i = 1; i < arraySize; i <<= 1)
            {
                for (int j = 0; j < i; j++)
                    reBitArray[j + i] = reBitArray[j] + arraySizeHarf;
                arraySizeHarf >>= 1;
            }
            return reBitArray;
        }


        #endregion

        public MainWindow()
        {
            InitializeComponent();

            yAxes.Maximum = 40;
            for (int i = 0; i < WaveIn.DeviceCount; i++)
            {
                var deviceInfo = WaveIn.GetCapabilities(i);
                Console.WriteLine(String.Format("Device {0}: {1}, {2} channels",
                    i, deviceInfo.ProductName, deviceInfo.Channels));
            }


            WaveIn waveIn = new WaveIn()
            {
                DeviceNumber = 0, // Default
            };
            waveIn.DataAvailable += WaveIn_DataAvailable;
            waveIn.WaveFormat = new WaveFormat(sampleRate: fs, channels: 1);
            waveIn.StartRecording();
        }

        private void WaveIn_DataAvailable(object sender, WaveInEventArgs e)
        {
            for (int index = 0; index < e.BytesRecorded; index += 2)
            {
                short sample = (short)((e.Buffer[index + 1] << 8) | e.Buffer[index + 0]);

                float sample32 = sample / (float)short.MaxValue;
                ProcessSample(sample32);
            }
        }

        List<double> _recorded = new List<double>(); // 音声データ
        public LineSeries _lineSeries = new LineSeries();




        private void ProcessSample(float sample)
        {
            var windowsize = size;
            _recorded.Add(sample);
            if (_recorded.Count == windowsize)
            {
                #region add
                // DFT用データ
                double[] dftIn = new double[size];
                double[] dftInIm = new double[size];
                DataPoint[] DftIn = new DataPoint[size];
                DataPoint[] DFTResult = new DataPoint[size];
                DataPoint[] FFTResult = new DataPoint[size];
                double[] data = new double[size];
                for (int i = 0; i < size; i++)
                {
                    dftInIm[i] = 0.0;
                }
                var window = MathNet.Numerics.Window.Hamming(windowsize);
                _recorded = _recorded.Select((v, i) => (double)v * window[i]).ToList();
                dftIn = _recorded.Take(size).ToArray();


                // FFT
                FFT(dftIn, dftInIm, out reFFT, out imFFT, (int)Math.Log(size, 2));
                // 波形显示
                for (int i = 0; i < size / 2; i++)
                {
                    if (i > 0)
                    {
                        float a = ((float)fs / size);
                        float x = (float)i * a;
                        double y = Math.Sqrt(reFFT[i] * reFFT[i] + imFFT[i] * imFFT[i]);
                        FFTResult[i] = new DataPoint(x, y);
                    }
                }

                line1.ItemsSource = FFTResult.Take((FFTResult.Count() / 20));
                #endregion
                _recorded.Clear();
            }
        }
    }
}

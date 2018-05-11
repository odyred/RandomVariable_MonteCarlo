﻿using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;
using System;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using Troschuetz.Random.Distributions.Continuous;

namespace ISAAR.MSolve.Analyzers
{
    public class SpectralRepresentation2DRandomField : IStochasticCoefficientsProvider
    {
        private readonly double XDirectionCorrelationLengthParameter;
        private readonly double YDirectionCorrelationLengthParameter;
        private readonly double spectrumStandardDeviation, cutoffError, frequencyIncrement, tolerance;
        private readonly int nPhi, frequencyIntervals, nptsVRF;
        private double kupper, period;
        private int frequencyCounter, frequencyCounterY;
        private double[] randomVariables = new double[0];
        private double[,] phi1,phi2;

        public SpectralRepresentation2DRandomField(double XDirectionCorrelationLengthParameter, double YDirectionCorrelationLengthParameter, double spectrumStandardDeviation, double cutoffError, 
            double frequencyIncrement = 0.1, int frequencyIntervals = 256, double tolerance = 1e-10)
        {
            this.XDirectionCorrelationLengthParameter = XDirectionCorrelationLengthParameter;
            this.YDirectionCorrelationLengthParameter = YDirectionCorrelationLengthParameter;
            this.spectrumStandardDeviation = spectrumStandardDeviation;
            this.cutoffError = cutoffError;
            this.frequencyIncrement = frequencyIncrement;
            this.frequencyIntervals = frequencyIntervals;
            this.tolerance = tolerance;
            Calculate();
        }

        public double SpectrumStandardDeviation { get { return spectrumStandardDeviation; } }
        public double Kupper { get { return kupper; } }
        public double Period { get { return period; } }
        public int NPtsVRF { get { return nptsVRF; } }
        public int FrequencyIntervals { get { return frequencyIntervals; } }
        public int NPhi { get { return nPhi; } }
        public int CurrentMCS { get; set; }
        public int CurrentFrequency { get; set; }

        private double AutoCorrelation(double tauX, double tauY)
        {
            return Math.Pow(spectrumStandardDeviation, 2) * Math.Exp(-Math.Abs(tauX) / XDirectionCorrelationLengthParameter - 
                Math.Abs(tauY) / YDirectionCorrelationLengthParameter);
        }

        private double SpectralDensity(double waveNumberX, double waveNumberY)
        {
            double bxby = XDirectionCorrelationLengthParameter * YDirectionCorrelationLengthParameter;
            double bxkx = waveNumberX * XDirectionCorrelationLengthParameter;
            double byky = waveNumberY * YDirectionCorrelationLengthParameter;
            return .25 / Math.PI * Math.Pow(spectrumStandardDeviation, 2) * bxby *
                Math.Exp(-.25 * (Math.Pow(bxkx, 2) + Math.Pow(byky, 2)));
        }

        private void Calculate()
        {
            double integral = AutoCorrelation(0,0) / 4d;
            double cumulativeSum = SpectralDensity(frequencyIncrement / 2d, frequencyIncrement / 2d) * Math.Pow(frequencyIncrement, 2);
            double rectangleVolumeX = 0;
            double rectangleVolumeY = 0;
            frequencyCounter = 2;

            while (cumulativeSum < (1d - cutoffError) * integral)
            {
                double ks = (frequencyCounter-1) * frequencyIncrement + frequencyIncrement / 2d;

                for (int i = 1; i < frequencyCounter; i++)
                {
                    double kx = i * frequencyIncrement + frequencyIncrement / 2d;
                    rectangleVolumeX += SpectralDensity(kx, ks) * Math.Pow(frequencyIncrement, 2);
                }

                for (int j = 1; j < frequencyCounter; j++)
                {
                    double ky = j * frequencyIncrement + frequencyIncrement / 2d;
                    rectangleVolumeY += SpectralDensity(ks, ky) * Math.Pow(frequencyIncrement, 2);
                }

                if ((rectangleVolumeX < tolerance) || (rectangleVolumeY<tolerance))
                    break;

                cumulativeSum += SpectralDensity(ks + frequencyIncrement/2d, ks + frequencyIncrement / 2d) * Math.Pow(frequencyIncrement, 2) + rectangleVolumeX + rectangleVolumeY;
                frequencyCounter++;
            }

            kupper = (frequencyCounter-1) * frequencyIncrement + frequencyIncrement/2d;
            double dk = kupper / (double)frequencyIntervals;
            period = 2d * Math.PI / dk;
        }

        public double GetCoefficient(double meanValue, double[] coordinate)
        {
            var dk = kupper / (double)frequencyIntervals;
            double randomCoefficient = 0;
            for (int i = 0; i < frequencyIntervals; i++)
            {
                for (int j = 0; j < frequencyIntervals; j++)
                {
                    randomCoefficient += Math.Sqrt(2) * (Math.Sqrt(2 * SpectralDensity(dk / 2 + i * dk, dk / 2 + j * dk) * dk) *
                                         (Math.Cos((dk / 2 + i * dk) * coordinate[0] + (dk / 2 + j * dk) * coordinate[1] + phi1[i,j])) +
                                          Math.Sqrt(2 * SpectralDensity(dk / 2 + i * dk, -(dk / 2 + j * dk)) * dk) *
                                          (Math.Cos((dk / 2 + i * dk) * coordinate[0] - (dk / 2 + j * dk) * coordinate[1] + phi2[i, j])));
                }

            }

            return randomCoefficient;

        }

        public double[] GetDerivative(double[] coordinate)
        {
            var dk = kupper / (double)frequencyIntervals;
            double[] derivative = {0, 0};
            for (int i = 0; i < frequencyIntervals; i++)
            {
                for (int j = 0; j < frequencyIntervals; j++)
                {
                    derivative[0] += -(dk / 2 + i * dk) * Math.Sqrt(2) * (Math.Sqrt(2 * SpectralDensity(dk / 2 + i * dk, dk / 2 + j * dk) * dk) *
                                                                          (Math.Sin((dk / 2 + i * dk) * coordinate[0] + (dk / 2 + j * dk) * coordinate[1] + phi1[i, j])) +
                                                                          Math.Sqrt(2 * SpectralDensity(dk / 2 + i * dk, -(dk / 2 + j * dk)) * dk) *
                                                                          (Math.Sin((dk / 2 + i * dk) * coordinate[0] - (dk / 2 + j * dk) * coordinate[1] + phi2[i, j])));
                    derivative[1] += -(dk / 2 + j * dk) * Math.Sqrt(2) * (Math.Sqrt(2 * SpectralDensity(dk / 2 + i * dk, dk / 2 + j * dk) * dk) *
                                                                          (Math.Sin((dk / 2 + i * dk) * coordinate[0] + (dk / 2 + j * dk) * coordinate[1] + phi1[i, j])) -
                                                                          Math.Sqrt(2 * SpectralDensity(dk / 2 + i * dk, -(dk / 2 + j * dk)) * dk) *
                                                                          (Math.Sin((dk / 2 + i * dk) * coordinate[0] - (dk / 2 + j * dk) * coordinate[1] + phi2[i, j])));
                }

            }

            return derivative;

        }

        public double[] RandomVariables
        {
            get
            {
                return randomVariables;
            }
            set
            {
                phi1 = new double[frequencyIntervals, frequencyIntervals];
                phi2 = new double[frequencyIntervals, frequencyIntervals];
                var Phi = new ContinuousUniformDistribution(0d, 2d * Math.PI);
                for (int i = 0; i < frequencyIntervals; i++)
                {
                    for (int j = 0; j < frequencyIntervals; j++)
                    {
                        phi1[i, j] = Phi.NextDouble();
                        phi2[i, j] = Phi.NextDouble();
                    }
                }
                randomVariables = value;
            }
        }
    }
}

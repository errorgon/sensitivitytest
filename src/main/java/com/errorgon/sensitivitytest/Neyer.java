package com.errorgon.sensitivitytest;

import com.errorgon.sensitivitytest.math.matrix.doubles.Matrix;
import com.errorgon.sensitivitytest.math.matrix.doubles.matrixtype.dense.DenseMatrix;
import com.errorgon.sensitivitytest.math.stats.regression.linear.glm.GLMProblem;
import com.errorgon.sensitivitytest.math.stats.regression.linear.glm.GeneralizedLinearModel;
import com.errorgon.sensitivitytest.math.stats.regression.linear.glm.distribution.Binomial;
import com.errorgon.sensitivitytest.math.stats.regression.linear.glm.distribution.link.Probit;
import com.errorgon.sensitivitytest.math.vector.doubles.dense.DenseVector;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Locale;

public class Neyer {

    String units;
    double muMin;
    double muMax;
    double sigmaGuess;

    double muMLE;
    double sigmaMLE;

    int failures;
    int successes;

    double minStimuli = Double.MAX_VALUE;
    double maxStimuli = Double.MIN_VALUE;

    double minSuccess = Double.MAX_VALUE;
    double maxFailure = Double.MIN_VALUE;

    double[] beta = new double[2];
    double betaMu;
    double betaSig;

    int precision;
    String neyerID;

    String response = "Linear normal response";

    ArrayList<Run> runList = new ArrayList<>();

    public Neyer(String units, double muMin, double muMax, double sigmaGuess) {
        this(units, muMin, muMax, sigmaGuess, 2);
    }

    public Neyer(String units, double muMin, double muMax, double sigmaGuess, int precision) {
        this(units, muMin, muMax, sigmaGuess, 2, new SimpleDateFormat().format(new SimpleDateFormat("yyyyMMMdd_HHmm", Locale.getDefault())));
    }

    public Neyer(String units, double muMin, double muMax, double sigmaGuess, int precision, String neyerID) {
        this.units = units;
        this.muMin = muMin;
        this.muMax = muMax;
        this.sigmaGuess = sigmaGuess;
        this.precision = precision;
        this.neyerID = neyerID;
    }


    public Run getRun() {
        // Phase 1
        if (runList.size() == 0) return firstRun();

        setMinMaxStimuli();

        int failures = 0;
        int successes = 0;
        for (int i = 0; i < runList.size(); i++) {
            if (!runList.get(i).getResult()) {
                failures++;
            } else if (runList.get(i).getResult()) {
                successes++;
            }
        }

        if (successes == 0) return allFailures();
        if (failures == 0) return allSuccesses();

        // Phase 2 -
        double diff = round(setMinMaxPassFail());
        if (diff > sigmaGuess) return diffSigma();

        // Compute MLE of Mu and Sigma
//        computeMLE();

        // If Sigma == 0, refine estimate of sigma
        if (sigmaMLE == 0) {
            return new Run(runList.size() + 1, (block4()), null);
        } else {
            return new Run(runList.size() + 1, (glmmle()), null);
        }
    }

    public void setRun(int trial, Double value, Boolean result) {
        if (result) successes++;
        else failures++;
        runList.add(new Run(trial, value, result));
    }

    public void setRun(Run run) {
        if (run.getResult()) successes++;
        else failures++;
        runList.add(run);
    }

    public ArrayList<Run> getRunList() { return runList; }
    public void setRunList(ArrayList<Run> runList) { this.runList = runList; }

    public int getFailures() { return failures; }
    public void setFailures(int failures) { this.failures = failures; }

    public int getSuccesses() { return successes; }
    public void setSuccesses(int successes) { this.successes = successes; }

    private Run firstRun() {
        return (new Run(runList.size() + 1, precision((muMax + muMin) / 2), null));
    }

    public String getUnits() {
        return units;
    }

    private void setMinMaxStimuli() {
        for (Run r : runList) {
            maxStimuli = Math.max(r.getValue(), maxStimuli);
            minStimuli = Math.min(r.getValue(), minStimuli);
        }
    }

    public String getNeyerID() {
        return neyerID;
    }

    public void setNeyerID(String neyerID) {
        this.neyerID = neyerID;
    }

    private double setMinMaxPassFail() {
        for (Run r : runList) {
            if (r.getResult()) minSuccess = Math.min(r.getValue(), minSuccess);
            else maxFailure = Math.max(r.getValue(), maxFailure);
        }
        return minSuccess - maxFailure;
    }

    private Run allFailures() {
        double first = precision((muMax + maxStimuli) / 2);   // 25
        double second = precision(maxStimuli - (2 * sigmaGuess));
        double third = precision(( 2 * maxStimuli) - minStimuli); //20

        return new Run(runList.size() + 1, Math.max(first, Math.max(second, third)), null);
    }

    private Run allSuccesses() {
        double first = precision((muMin + minStimuli) / 2);
        double second = precision(minStimuli - (2 * sigmaGuess));
        double third = precision((2 * minStimuli) - maxStimuli);

        return new Run(runList.size() + 1, Math.min(first, Math.min(second, third)), null);
    }

    private Run diffSigma() {
        return new Run(runList.size() + 1, precision((maxFailure + minSuccess) / 2), null);
    }

    private void computeMLE() {
        int n = 0;
        double sum = 0;
        for (Run r : runList) {
            if (r.getResult()) {
                n++;
                sum += r.getValue();
            }
        }
        muMLE = (sum / n);

        double add = 0;
        for (Run r : runList) {
            if (r.getResult()) {
                add += Math.pow((r.getValue() - muMLE), 2);
            }
        }
        sigmaMLE = Math.sqrt(add / n);
    }

    private double[][] getFisher(double sigmaEstimate) {
        double fisher[][] = new double[2][2];

//        n = count
        int count = runList.size();
        for (Run r : runList) {
            double k = 0;
            if (r.getResult()) {
                count++;
                k = (r.getValue() - muMLE) / sigmaMLE;
                double p = pnorm(k) * (1 - pnorm(k));
                double z = dnorm(k);
                double v = (count * Math.pow(z, 2)) / p;
            }
        }

        return fisher;
    }

    public double Sk(double k, double[][] b) {
        return b[0][0] * Math.pow(k, 2) - 2.0 * b[0][1] * k + b[1][1];
    }

    public double Gk(double k) {
        double pk = pnorm(k);
        double gk = dnorm(k) / Math.sqrt(pk * (1 - pk));
        return Math.pow(gk, 2);
    }

    // Equiv to R pnorm
    // return the integral from -inf to q of the pdf of the normal distribution
    // where q is a z-score
    // returns the cumulative distribution function (CDF)
    public double pnorm(double q) {
        NormalDistribution nd = new NormalDistribution();
        return nd.cumulativeProbability(q);
    }

    // Equiv to R dnorm
    // returns the value of the probability density function for the
    // normal distribution
    public double dnorm(double value) {
        NormalDistribution nd = new NormalDistribution();
        return nd.density(value);
    }

    // derivative of g^2*s, where g=dk/sqrt(pk*(1-pk))
    // s=b11*k^2-2*b12*k+b22, b=Information Matrix
    public double derivativeGs(double k, double [][] matrix) {
        double pk = pnorm(k);
        double dk = dnorm(k);
        double sk = matrix[0][0] * Math.pow(k, 2) - 2 * matrix[0][1] * k + matrix[1][1];
        double j = 2 * (matrix[0][0] - sk) * k - (2 * matrix[0][1] + dk * sk * (1 - 2 * pk) / (pk * (1 - pk)));
        return j;
    }

    // b=information matrix; presumption: kmax & b12 have opposite signs
    // first part finds [l2,0] or [0,l2] containing the max g(k)^2*s(k)
    // second part solves for the zero of d/dk of g(k)^2*s(k)
    // this function uses functions Sk, Gk and dgs
    public double kstar(double[][] b) {
        double del = -1.0;
        double k1 = 0.0;
        double val1 = Gk(k1) * Sk(k1, b);
        double val2 = val1 + 1;

        if (b[0][1] <= 0) {
            del = 1;
        }

        double k2 = 0.0;
        while (val2 > val1) {
            k2 = k1 + del;
            val2 = Gk(k2) * Sk(k2, b);
            k1 = k2;
        }

        double eps = .000001;
        k1 = 0;
        double v1 = derivativeGs(k1, b);
        double v2 = derivativeGs(k2, b);

        while (Math.abs(k2 - k1) > eps | Math.abs(v2 - v1) > eps) {
            double kmid = (k1 + k2) / 2;
            double vmid = (derivativeGs(kmid, b));
            if (v1 * vmid > 0) {
                v1 = vmid;
                k1 = kmid;
            } else {
                v2 = vmid;
                k2 = kmid;
            }
        }
        double kmax = (k1 + k2) / 2;
        return kmax;
    }

    // block == 4 equals sigma0
    private double block4() {
//        m = (m1 + M0) / 2;
        muMLE = round((minSuccess + maxFailure) / 2);

//        es = sg;
        sigmaMLE = sigmaGuess;

//        sg = .8 * sg;
        sigmaGuess = round(0.8 * sigmaGuess);

//        m = max(xlo, min(m, xhi));
        muMLE = Math.max(minStimuli, Math.min(muMLE, maxStimuli));

//        es = min(es, (xhi - xlo));
        sigmaMLE = Math.min(sigmaMLE, maxStimuli - minStimuli);

//        b = yinfomat(d0, m, es)$infm ;
        double[][] matrix = yinformat();

        double xNext = muMLE + kstar(matrix) * sigmaMLE;

        return precision(xNext);
    }

    public double[][] yinformat(double mu, double sigma) {
        int count = runList.size();

        for (Run r : runList) {
            r.k = (r.getValue() - mu) /sigma;
            r.p = pnorm(r.k) * (1 - pnorm(r.k));
            r.z = dnorm(r.k);
            r.v = Math.pow(r.z, 2) / r.p;
        }

        double sumV = 0;
        for (Run r : runList) {
            sumV += r.v;
        }

        double sumVK = 0;
        for (Run r : runList) {
            sumVK += (r.v * r.k);
        }

        double sumVK2 = 0;
        for (Run r : runList) {
            sumVK2 += (r.v * Math.pow(r.k, 2));
        }

        double[][] informationMatrix = new double[2][2];
        informationMatrix[0][0] = sumV;
        informationMatrix[1][0] = informationMatrix[0][1] = sumVK;
        informationMatrix[1][1] = sumVK2;

        double determinant = determinant(informationMatrix);

        double[][] inverse = inverseMatrix(informationMatrix);

        double correlationCoefficient = correlationCoefficient(inverse);

        return informationMatrix;
    }

    public double[][] yinformat() {
        return yinformat(muMLE, sigmaMLE);
    }

    private double glmmle() {

        double[][] input = new double[1][runList.size()];
        double[] result = new double[runList.size()];

        for (int i = 0; i < runList.size(); i++) {
            input[0][i] = runList.get(i).getValue();
            if (runList.get(i).getResult()) {
                result[i] = 1;
            } else {
                result[i] = 0;
            }
        }

        Matrix Xt = new DenseMatrix(input);
        GLMProblem problem = new GLMProblem(
                new DenseVector(result),
                Xt.t(),
                true, new Binomial(new Probit()));

        GeneralizedLinearModel instance = new GeneralizedLinearModel(problem);

        beta[0] = instance.beta.betaHat.toArray()[0];
        beta[1] = instance.beta.betaHat.toArray()[1];

        double mu = -instance.beta.betaHat.toArray()[1] / instance.beta.betaHat.toArray()[0];
        betaMu = mu;
        double sigma = 1 / instance.beta.betaHat.toArray()[0];
        betaSig = sigma;
        double minX = Double.MAX_VALUE;
        double maxX = Double.MIN_VALUE;
        for (int i = 0; i < input[0].length; i++) {
            if (input[0][i] < minX) {
                minX = input[0][i];
            }
            if (input[0][i] > maxX) {
                maxX = input[0][i];
            }
        }

        double mu4 = Math.max(minX, Math.min(mu, maxX));

        double sg4 = Math.min(sigma, maxX - minX);

        double[][] b = yinformat(mu4, sg4);

        double xNext = mu4 + kstar(b) * sg4;

        return precision(xNext);
    }

    private double determinant(double[][] matrix) {
        return (matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]);
    }

    private double[][] inverseMatrix(double[][] matrix) {
        double det = 1 / determinant(matrix);
        double[][] inverse = new double[2][2];
        inverse[0][0] = matrix[1][1] * det;
        inverse[1][1] = matrix[0][0] * det;
        inverse[0][1] = -matrix[0][1] * det;
        inverse[1][0] = -matrix[1][0] * det;
        return inverse;
    }

    private double correlationCoefficient(double[][] matrix) {
        return matrix[0][1] / (Math.sqrt(matrix[0][0] * matrix[1][1]));
    }

    double round(double num) {
        double r = Math.pow(100000, 2);
        return Math.round(num * r) / r;
    }

    double precision(double num) {
        double r = Math.pow(10, precision);
        return Math.round(num * r) / r;
    }

    public double[] getBeta() {
        return beta;
    }

    public double getBetaMu() {
        return betaMu;
    }

    public double getBetaSig() {
        return betaSig;
    }

    public class Run {
        int trial;
        Double value;
        Boolean result;

        double k = 0;
        double p = 0;
        double z = 0;
        double v = 0;


        public Run(int trial, Double value, Boolean result) {

            this.trial = trial;
            this.value = value;
            this.result = result;
        }

        public String getUnits() { return Neyer.this.units; }
        public int getPrecision() { return Neyer.this.precision;}

        public int getTrial() { return trial; }
        public void setTrial(int trial) { this.trial = trial; }

        public Double getValue() { return value; }
        public void setValue(Double value) { this.value = value; }

        public Boolean getResult() { return result; }
        public void setResult(Boolean result) { this.result = result; }

        @Override
        public String toString() {
            return "Trial# " + trial + " Value: " + value + " Result: " + result;
        }
    }

}

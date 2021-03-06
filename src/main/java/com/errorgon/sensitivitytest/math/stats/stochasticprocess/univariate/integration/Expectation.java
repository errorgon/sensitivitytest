/*
 * Copyright (c) Numerical Method Inc.
 * http://www.numericalmethod.com/
 * 
 * THIS SOFTWARE IS LICENSED, NOT SOLD.
 * 
 * YOU MAY USE THIS SOFTWARE ONLY AS DESCRIBED IN THE LICENSE.
 * IF YOU ARE NOT AWARE OF AND/OR DO NOT AGREE TO THE TERMS OF THE LICENSE,
 * DO NOT USE THIS SOFTWARE.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITH NO WARRANTY WHATSOEVER,
 * EITHER EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION,
 * ANY WARRANTIES OF ACCURACY, ACCESSIBILITY, COMPLETENESS,
 * FITNESS FOR A PARTICULAR PURPOSE, MERCHANTABILITY, NON-INFRINGEMENT, 
 * TITLE AND USEFULNESS.
 * 
 * IN NO EVENT AND UNDER NO LEGAL THEORY,
 * WHETHER IN ACTION, CONTRACT, NEGLIGENCE, TORT, OR OTHERWISE,
 * SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 * ANY CLAIMS, DAMAGES OR OTHER LIABILITIES,
 * ARISING AS A RESULT OF USING OR OTHER DEALINGS IN THE SOFTWARE.
 */
package com.errorgon.sensitivitytest.math.stats.stochasticprocess.univariate.integration;

import com.errorgon.sensitivitytest.math.stats.descriptive.moment.Mean;
import com.errorgon.sensitivitytest.math.stats.descriptive.moment.Variance;
import com.errorgon.sensitivitytest.math.stats.stochasticprocess.timepoints.EvenlySpacedGrid;
import com.errorgon.sensitivitytest.math.stats.stochasticprocess.timepoints.TimeGrid;
import com.errorgon.sensitivitytest.math.stats.stochasticprocess.univariate.brownian.RandomWalk;
import com.errorgon.sensitivitytest.math.stats.stochasticprocess.univariate.integration.sde.Construction;

/**
 * This class computes the expectation of the following class of integrals.
 *
 * <blockquote><code><pre>
 *    /
 * E( | f(B(t)) dt )
 *    /
 * </pre></code></blockquote>
 *
 * <p>
 * <i>f</i> is not necessarily an adapted function.
 *
 * <p>
 * The Brownian paths are generated by Monte Carlo simulations.
 *
 * @author Haksun Li
 */
public class Expectation {

    /**
     * the integrator to compute the integral for each filtration/path
     */
    public final Integrator I;
    /**
     * the beginning time of the integral time interval
     */
    public final double t0;
    /**
     * the ending time of the integral time interval
     */
    public final double t1;
    /**
     * the number of discretization in the integral time interval
     */
    public final int n;
    /**
     * the number of simulations
     */
    public final int nSim;
    /**
     * the mean of the stochastic integral
     */
    private Mean mean = new Mean();
    /**
     * the variance of the stochastic integral
     */
    private Variance var = new Variance();

    /**
     * Compute the expectation for the integral of a stochastic process.
     * 
     * @param I the integral of a stochastic process
     * @param t0 the beginning time of the integral time interval
     * @param t1 the ending time of the integral time interval
     * @param n the number of discretization in the integral time interval
     * @param nSim the number of simulations
     */
    public Expectation(Integrator I, double t0, double t1, int n, int nSim) {
        this.I = I;
        this.t0 = t0;
        this.t1 = t1;
        this.n = n;
        this.nSim = nSim;

        simulate();
    }

    /**
     * Compute the mean of the integral.
     *
     * @return the mean of the integral
     */
    public double mean() {
        return mean.value();
    }

    /**
     * Compute the variance of the integral.
     *
     * @return the variance of the integral
     */
    public double var() {
        return var.value();
    }

    /**
     * Compute the integral values for a set of simulated filtration/paths.
     * Then we take the mean and variance.
     */
    private void simulate() {
        TimeGrid T = new EvenlySpacedGrid(t0, t1, n);

        for (int i = 0; i < nSim; ++i) {
            Construction B = new RandomWalk(T);
            Filtration FT = new Filtration(B.nextRealization(0));

            double value = I.integral(FT);

            mean.addData(value);
            var.addData(value);
        }
    }
}

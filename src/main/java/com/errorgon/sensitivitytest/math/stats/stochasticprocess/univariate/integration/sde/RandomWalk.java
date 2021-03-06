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
package com.errorgon.sensitivitytest.math.stats.stochasticprocess.univariate.integration.sde;

import com.errorgon.sensitivitytest.math.stats.random.univariate.RandomLongGenerator;
import com.errorgon.sensitivitytest.math.stats.random.univariate.uniform.MersenneTwister;
import com.errorgon.sensitivitytest.math.stats.stochasticprocess.timepoints.TimeGrid;
import com.errorgon.sensitivitytest.math.stats.stochasticprocess.univariate.sde.DiscretizedSDE;
import com.errorgon.sensitivitytest.math.stats.stochasticprocess.univariate.sde.Ft;

/**
 * This is the Random Walk construction of a stochastic process per SDE specification.
 *
 * <p>
 * Given an SDE, we compute the increments of the stochastic process.
 * A realization is the cumulative summation of these independent increments.
 * 
 * @author Haksun Li
 */
public class RandomWalk implements Construction {

    /**
     * the SDE specification, in discretized form
     */
    public final DiscretizedSDE sde;
    /**
     * the set of discretized time points
     */
    public final TimeGrid timePoints;
    /**
     * generate the seeds for the iterators
     */
    private RandomLongGenerator uniform = new MersenneTwister();//TODO: make this a pass-in argument

    public class Realization implements com.errorgon.sensitivitytest.math.stats.stochasticprocess.univariate.Realization {

        /**
         * the ID of this particular realization
         */
        public final long id;
        /**
         * the initial value of the realization
         */
        public final double x0;

        private Realization(double x0) {
            this.id = uniform.nextLong();
            this.x0 = x0;
        }

        public int size() {
            return timePoints.size();
        }

        public Iterator iterator() {
            return new Realization.Iterator(size(), id) {

                private double t0 = 0, t1;
                private double xt = x0;// initialize the starting value
                private Ft ft = sde.getNewFt();

                public double t(int index) {
                    return timePoints.t(index);
                }

                /**
                 * {@inheritDoc}
                 *
                 * <p>
                 * This is an implementation of the Random Walk construction.
                 *
                 * @param index the index to the current <tt>Entry</tt>
                 * @return the current value
                 */
                @Override
                public double xt(int index) {
                    t1 = t(index);

                    ft.setDt(t1 - t0);
                    ft.setXt(xt);
                    ft.setZt(Zt());

                    double dx = sde.dXt(ft);
                    xt += dx;
                    t0 = t1;

                    return xt;
                }

                public void remove() {
                    throw new UnsupportedOperationException("time series is immutable");
                }
            };
        }

        public double[] toArray() {
            Realization.Iterator it = iterator();

            double[] wt = new double[size()];

            for (int i = 0; it.hasNext(); ++i) {
                wt[i] = it.nextValue();
            }

            return wt;
        }

        public double lastValue() {
            Realization.Iterator it = iterator();

            double wt = Double.NaN;

            for (int i = 0; it.hasNext(); ++i) {
                wt = it.nextValue();
            }

            return wt;
        }
    }

    /**
     * Construct a univariate stochastic process from an SDE.
     * The realizations are generated by the Random Walk method.
     *
     * @param sde an SDE
     * @param timePoints specifying the time points in a grid
     */
    public RandomWalk(DiscretizedSDE sde, TimeGrid timePoints) {
        this.sde = sde;
        this.timePoints = timePoints;
    }

    public Realization nextRealization(double x0) {
        return new Realization(x0);
    }

    public void seed(long seed) {
        uniform.seed(seed);
    }
}

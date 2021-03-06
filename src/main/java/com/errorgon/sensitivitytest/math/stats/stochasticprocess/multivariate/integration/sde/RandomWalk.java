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
package com.errorgon.sensitivitytest.math.stats.stochasticprocess.multivariate.integration.sde;

import com.errorgon.sensitivitytest.math.stats.random.univariate.RandomLongGenerator;
import com.errorgon.sensitivitytest.math.stats.random.univariate.uniform.MersenneTwister;
import com.errorgon.sensitivitytest.math.stats.stochasticprocess.multivariate.sde.Ft;
import com.errorgon.sensitivitytest.math.matrix.doubles.Matrix;
import com.errorgon.sensitivitytest.math.stats.stochasticprocess.multivariate.sde.DiscretizedSDE;
import com.errorgon.sensitivitytest.math.stats.stochasticprocess.timepoints.TimeGrid;
import com.errorgon.sensitivitytest.math.vector.doubles.Vector;
import static com.errorgon.sensitivitytest.math.matrix.doubles.operation.CreateMatrix.rbind;

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

    /*
     * the realization, constructed by the Random Walk method, of a multivariate stochastic process
     */
    public class MultiVariateRealization implements com.errorgon.sensitivitytest.math.stats.stochasticprocess.multivariate.MultiVariateRealization {

        /**
         * the ID of this particular realization
         */
        public final long id;
        /**
         * the initial value of the realization
         */
        public final Vector x0;

        private MultiVariateRealization(Vector x0) {
            this.id = uniform.nextLong();
            this.x0 = x0;
        }

        public int size() {
            return timePoints.size();
        }

        public int dimension() {
            return x0.size();
        }

        public Iterator iterator() {
            return new MultiVariateRealization.Iterator(sde.nB(), size(), id) {

                private double t0 = 0, t1;
                private Vector xt = x0;// initialize the starting value
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
                public Vector xt(int index) {
                    t1 = t(index);

                    ft.setDt(t1 - t0);
                    ft.setXt(xt);
                    ft.setZt(Zt());

                    Vector dx = sde.dXt(ft);
                    xt = xt.add(dx);
                    t0 = t1;

                    return xt;
                }

                public void remove() {
                    throw new UnsupportedOperationException("Not supported yet.");
                }
            };
        }

        public Vector lastValue() {
            MultiVariateRealization.Iterator it = iterator();

            Vector wt = null;

            for (int i = 0; it.hasNext(); ++i) {
                wt = it.nextValue();
            }

            return wt;
        }

        public Matrix toMatrix() {
            MultiVariateRealization.Iterator it = iterator();

            Vector[] wt = new Vector[size()];

            for (int i = 0; it.hasNext(); ++i) {
                wt[i] = it.nextValue();
            }

            return rbind(wt);
        }
    }

    /**
     * Construct a multivariate stochastic process from an SDE.
     * The realizations are generated by the Random Walk method.
     *
     * @param sde an SDE
     * @param timePoints specifying the time points in a grid
     */
    public RandomWalk(DiscretizedSDE sde, TimeGrid timePoints) {
        this.sde = sde;
        this.timePoints = timePoints;
    }

    public MultiVariateRealization nextRealization(Vector x0) {
        return new MultiVariateRealization(x0);
    }

    public void seed(long seed) {
        uniform.seed(seed);
    }
}

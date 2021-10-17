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
package com.errorgon.sensitivitytest.math.optimization.geneticalgorithm.minimizer.local;

import com.errorgon.sensitivitytest.math.analysis.function.rn2r1.RealScalarFunction;
import com.errorgon.sensitivitytest.math.optimization.Minimizer;
import com.errorgon.sensitivitytest.math.optimization.geneticalgorithm.minimizer.simplegrid.SimpleCellFactory;
import com.errorgon.sensitivitytest.math.optimization.problem.C2OptimProblemImpl;
import com.errorgon.sensitivitytest.math.optimization.problem.IterativeMinimizer;
import com.errorgon.sensitivitytest.math.optimization.problem.OptimProblem;
import com.errorgon.sensitivitytest.math.stats.random.univariate.RandomLongGenerator;
import com.errorgon.sensitivitytest.math.vector.doubles.Vector;

/**
 * A {@code LocalSearchCellFactory} produces {@code LocalSearchCellFactory.LocalSearchCell}s.
 * A {@code LocalSearchCellFactory.LocalSearchCell} is a chromosome for a real valued function (an optimization problem) and a candidate solution.
 * Its mutation operation is by a local minimization method, Nelder-Mead, BFGS, searching the neighborhood of the current solution.
 *
 * @param <P> the optimization problem type
 * @param <T> the minimizer type for local search
 * @author Haksun Li
 */
public class LocalSearchCellFactory<P extends OptimProblem, T extends Minimizer<OptimProblem, IterativeMinimizer<Vector>>> extends SimpleCellFactory {

    /**
     * This factory constructs a new {@code Minimizer} for each mutation operation.
     *
     * @param <U> the minimizer type for local search
     */
    public static interface MinimizerFactory<U extends Minimizer<OptimProblem, IterativeMinimizer<Vector>>> {

        /**
         * Construct a new instance of {@code Minimizer} for a mutation operation.
         *
         * @return a new instance of {@code Minimizer}
         */
        public U newInstance();
    }

    private final MinimizerFactory<T> factory;

    /**
     * Construct an instance of a {@code LocalSearchCellFactory}.
     *
     * @param factory a factory that generates a new instance of a {@link Minimizer} for local search
     * @param uniform a uniform random number generator
     */
    public LocalSearchCellFactory(MinimizerFactory<T> factory, RandomLongGenerator uniform) {
        super(0.1, uniform);
        this.factory = factory;
    }

    @Override
    public LocalSearchCell getSimpleCell(RealScalarFunction f, Vector x) {
        return new LocalSearchCell(f, x);
    }

    /**
     * A {@code LocalSearchCell} implements the two genetic operations.
     * <ul>
     * <li>Mutation by a local search in the neighborhood;
     * <li>Crossover by taking the midpoint (average) of two cells.
     * </ul>
     */
    public class LocalSearchCell extends SimpleCellFactory.SimpleCell {

        protected LocalSearchCell(RealScalarFunction f, Vector x) {
            super(f, x);
        }

        /**
         * Mutate by a local search in the neighborhood.
         *
         * @return a mutant chromosome
         */
        @Override
        public LocalSearchCell mutate() {
            try {
                T solver = factory.newInstance();//need a new minimizer in each mutation to do parallelization
                IterativeMinimizer<Vector> soln = solver.solve(new C2OptimProblemImpl(f()));
                Vector z = soln.search(x());

                return getSimpleCell(f(), z);
            } catch (Exception ex) {
                return getSimpleCell(f(), x());//unchanged
            }
        }
    }

}

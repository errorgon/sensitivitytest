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
package com.errorgon.sensitivitytest.math.optimization.geneticalgorithm.minimizer.simplegrid;

import com.errorgon.sensitivitytest.math.Constant;
import com.errorgon.sensitivitytest.math.analysis.function.rn2r1.RealScalarFunction;
import com.errorgon.sensitivitytest.math.interval.RealInterval;
import com.errorgon.sensitivitytest.math.optimization.Minimizer;
import com.errorgon.sensitivitytest.math.optimization.geneticalgorithm.Chromosome;
import com.errorgon.sensitivitytest.math.optimization.geneticalgorithm.GeneticAlgorithm;
import com.errorgon.sensitivitytest.math.optimization.initialization.UniformDistributionOverBox2;
import com.errorgon.sensitivitytest.math.optimization.problem.IterativeMinimizer;
import com.errorgon.sensitivitytest.math.optimization.problem.OptimProblem;
import com.errorgon.sensitivitytest.math.stats.random.univariate.RandomLongGenerator;
import com.errorgon.sensitivitytest.math.stats.random.univariate.uniform.UniformRng;
import com.errorgon.sensitivitytest.math.vector.doubles.ImmutableVector;
import com.errorgon.sensitivitytest.math.vector.doubles.Vector;
import static java.lang.Math.*;
import java.util.ArrayList;

/**
 * This minimizer is a simple global optimization method.
 * It puts a mesh over the feasible region and then locally searches (optimizes) the neighborhood around each mesh point.
 * The algorithm tries to escape the local minimums by crossing over other local minimums using a genetic algorithm.
 *
 * @author Haksun Li
 */
public class SimpleGridMinimizer implements Minimizer<OptimProblem, IterativeMinimizer<Vector>> {

    /**
     * This factory constructs a new {@code SimpleCellFactory} for each minimization problem.
     */
    public static interface NewCellFactoryCtor {

        /**
         * Construct a new instance of {@code SimpleCellFactory} for a minimization problem.
         *
         * @return a new instance of {@code SimpleCellFactory}
         */
        public SimpleCellFactory newCellFactory();
    }

    /**
     * This is the solution to a minimization problem using {@code SimpleGridMinimizer}.
     */
    protected class Solution extends GeneticAlgorithm implements IterativeMinimizer<Vector> {

        protected Vector[] initials;
        protected int iteration = 0;
        protected int nNoChanges = 0;
        protected double fminLast = Double.POSITIVE_INFINITY;
        protected double fmin = Double.POSITIVE_INFINITY;
        protected Vector xmin;
        protected final RealScalarFunction f;
        protected final SimpleCellFactory factory;

        protected Solution(RealScalarFunction f) {
            super(SimpleGridMinimizer.this.parallel, SimpleGridMinimizer.this.uniform);
            this.f = f;
            this.factory = factoryCtor.newCellFactory();
        }

        /**
         * The initial population is generated by putting a uniform mesh/grid/net over the entire region.
         * The grid points are the chromosomes in the first population.
         * The population size is proportional to the number of available cores.
         * The region bounds are determined from the initial guesses.
         *
         * @return the first population
         */
        @Override
        protected ArrayList<Chromosome> initialization() {
            final int dim = f.dimensionOfDomain();

            //handle the case there is only initial when bounds cannot be determined
            Vector[] initials = this.initials;
            if (initials.length == 1) {
                initials = new Vector[3];
                initials[0] = this.initials[0];
                initials[1] = this.initials[0].scaled(0.5);
                initials[2] = this.initials[0].scaled(1.5);
            }

            //find bounds
            RealInterval[] bounds = new RealInterval[dim];
            for (int i = 0; i < dim; ++i) {
                double lower = Double.POSITIVE_INFINITY;
                double upper = Double.NEGATIVE_INFINITY;
                for (int j = 0; j < initials.length; ++j) {
                    if (initials[j].get(i + 1) < lower) {
                        lower = initials[j].get(i + 1);
                    }

                    if (initials[j].get(i + 1) > upper) {
                        upper = initials[j].get(i + 1);
                    }
                }

                bounds[i] = new RealInterval(lower, upper);
            }

            //getInitials a more or less uniform mesh over the region
            int nProcessors = Runtime.getRuntime().availableProcessors();//(d-1)^dim = #processors
            int discretization = (int) exp(log(nProcessors) / dim) + 1;
            discretization = Math.max(discretization, 4);
            UniformDistributionOverBox2 buildInitials = new UniformDistributionOverBox2(0.1, bounds, discretization, uniform);

            Vector[] mesh = buildInitials.getInitials();
            ArrayList<Chromosome> cells = new ArrayList<Chromosome>(mesh.length);
            for (int i = 0; i < mesh.length; ++i) {
                cells.add(factory.getSimpleCell(f, mesh[i]));
            }

            return cells;
        }

        /**
         * This genetic algorithm terminates if
         * <ul>
         * <li>the minimum does not improve for a fixed number of iterations, or
         * <li>the maximum number of iterations is exceeded.
         * </ul>
         *
         * @return {@code true} if a convergence is found
         */
        @Override
        protected boolean isConverged() {
            double min = minimum();
            double dmin = abs(min - fminLast);
            if (min < fminLast) {
                fminLast = min;
            }

            if (dmin < epsilon * abs(min) + Constant.EPSILON) {
                nNoChanges++;
            } else {
                nNoChanges = 0;
            }

            if (nNoChanges > nStableIterations) {
                return true;
            }

            if (iteration >= maxIterations) {
                return true;
            }

            return false;
        }

        @Override
        public void setInitials(Vector... initials) {
            this.initials = initials;
        }

        @Override
        public Object step() {
            ++iteration;
            super.step();

            //update state info
            SimpleCellFactory.SimpleCell best = (SimpleCellFactory.SimpleCell) getBest(0);
            xmin = best.x();
            fmin = getBest(0).fitness();

            return true;
        }

        @Override
        public Vector search(Vector... initials) {
            setInitials(initials);
            super.run();

            return minimizer();
        }

        @Override
        public double minimum() {
            return fmin;
        }

        @Override
        public ImmutableVector minimizer() {
            return new ImmutableVector(xmin);
        }
    }

    protected final boolean parallel;
    protected final RandomLongGenerator uniform;
    protected final NewCellFactoryCtor factoryCtor;
    protected final double epsilon;
    protected final int maxIterations;
    protected final int nStableIterations;

    /**
     * Construct a {@code SimpleGridMinimizer} to solve unconstrained minimization problems.
     *
     * @param factoryCtor       a factory that constructs a new instance of {@code SimpleCellFactory} for each problem
     * @param parallel          {@code true} if the algorithm is to run in parallel (multi-core)
     * @param uniform           a uniform random number generator
     * @param epsilon           a precision parameter: when a number |x| ≤ ε, it is considered 0
     * @param maxIterations     the maximum number of iterations
     * @param nStableIterations The solution is considered converged if the minimum does not change over this many iterations.
     */
    public SimpleGridMinimizer(NewCellFactoryCtor factoryCtor, boolean parallel, RandomLongGenerator uniform, double epsilon, int maxIterations, int nStableIterations) {
        this.parallel = parallel;
        this.uniform = uniform;
        this.factoryCtor = factoryCtor;
        this.epsilon = epsilon;
        this.maxIterations = maxIterations;
        this.nStableIterations = nStableIterations;
    }

    /**
     * Construct a {@code SimpleGridMinimizer} to solve unconstrained minimization problems.
     *
     * @param parallel      {@code true} if the algorithm is to run in parallel (multi-core)
     * @param epsilon       a precision parameter: when a number |x| ≤ ε, it is considered 0
     * @param maxIterations the maximum number of iterations
     */
    public SimpleGridMinimizer(boolean parallel, double epsilon, int maxIterations) {
        this(
                new NewCellFactoryCtor() {

                    @Override
                    public SimpleCellFactory newCellFactory() {
                        return new SimpleCellFactory(0.1, new UniformRng());
                    }
                },
                parallel, new UniformRng(), epsilon, maxIterations, 10);
    }

    @Override
    public IterativeMinimizer<Vector> solve(OptimProblem problem) throws Exception {
        return new Solution(problem.f());
    }
}

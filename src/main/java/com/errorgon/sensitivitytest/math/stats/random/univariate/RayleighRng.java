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
package com.errorgon.sensitivitytest.math.stats.random.univariate;

import com.errorgon.sensitivitytest.math.stats.random.univariate.uniform.UniformRng;

/**
 * This random number generator samples from the Rayleigh distribution using the inverse transform sampling method.
 *
 * @author Haksun Li
 * @see <a href="http://en.wikipedia.org/wiki/Rayleigh_distribution">Wikipedia: RayleighDistribution distribution</a>
 */
public class RayleighRng extends InverseTransformSampling {

    /**
     * Construct a random number generator to sample from the Rayleigh distribution.
     *
     * @param sigma the standard deviation
     * @param rng   a uniform random number generator
     */
    public RayleighRng(double sigma, RandomLongGenerator rng) {
        super(new com.errorgon.sensitivitytest.math.stats.distribution.univariate.RayleighDistribution(sigma), rng);
    }

    /**
     * Construct a random number generator to sample from the Rayleigh distribution.
     *
     * @param sigma the standard deviation
     */
    public RayleighRng(double sigma) {
        super(new com.errorgon.sensitivitytest.math.stats.distribution.univariate.RayleighDistribution(sigma), new UniformRng());
    }
}

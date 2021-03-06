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
package com.errorgon.sensitivitytest.math.stats.stochasticprocess.multivariate.sde.coefficients;

import com.errorgon.sensitivitytest.math.stats.stochasticprocess.multivariate.sde.Ft;
import com.errorgon.sensitivitytest.math.matrix.doubles.Matrix;

/**
 * This represents the diffusion term, <i>σ</i>, of an SDE. It is of this form: <i>σ(dt, Xt, Zt, ...)</i>.
 *
 * <p>
 * Note that we are passing in the time differential, <i>dt</i>, instead of the time itself.
 * If we want to compute for a time <i>t</i>, the subclass would need to accumulate the <i>dt</i>s explicitly.
 * 
 * @author Haksun Li
 *
 * @see "Fima C. Klebaner. Introduction to Stochastic Calculus with Applications. 2nd ed. Section 4.7. Imperial College Press. 2006."
 */
public interface Diffusion {

    /**
     * <i>σ(dt, Xt, Zt, ...)</i>
     * 
     * @param ft filtration
     * @return a diffusion matrix
     */
    public Matrix evaluate(Ft ft);

    /**
     * Get the dimension of the process.
     *
     * @return the number of rows in the diffusion matrix
     */
    public int nRows();

    /**
     * Get the number of independent Brownian motions.
     * 
     * @return the number of columns in the diffusion matrix
     */
    public int nCols();
}

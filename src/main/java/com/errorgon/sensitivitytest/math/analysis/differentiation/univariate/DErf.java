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
package com.errorgon.sensitivitytest.math.analysis.differentiation.univariate;

import static com.errorgon.sensitivitytest.math.Constant.ROOT_PI;
import com.errorgon.sensitivitytest.math.analysis.function.rn2r1.univariate.UnivariateRealFunction;
import com.errorgon.sensitivitytest.math.analysis.function.special.gaussian.Erf;
import static java.lang.Math.exp;

/**
 * This is the first order derivative function of the Error function, {@link Erf}.
 * \[
 * {d \over dx}\mbox{Erf}(x) = \frac{2}{\sqrt{\pi}}\exp(-x^2)
 * \]
 *
 * @author Haksun Li
 * @see Erf
 */
public class DErf extends UnivariateRealFunction {

    @Override
    public double evaluate(double x) {
        return 2d / ROOT_PI * exp(-x * x);
    }
}

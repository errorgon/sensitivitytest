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
package com.errorgon.sensitivitytest.math.analysis.function.rn2r1.univariate;

import com.errorgon.sensitivitytest.math.analysis.function.rn2r1.RealScalarFunction;
import com.errorgon.sensitivitytest.math.vector.doubles.Vector;

/**
 * A univariate real function takes one real argument and outputs one real value.
 * That is, <i>y = f(x)</i>.
 *
 * @author Haksun Li
 */
public abstract class UnivariateRealFunction implements RealScalarFunction {

    @Override
    public int dimensionOfDomain() {
        return 1;
    }

    @Override
    public int dimensionOfRange() {
        return 1;
    }

    @Override
    public Double evaluate(Vector x) {
        if (x.size() != 1) {
            throw new EvaluationException("this is a univariate function");
        }
        return evaluate(x.get(1));
    }

    /**
     * Evaluate <i>y = f(x)</i>.
     *
     * @param x <i>x</i>
     * @return <i>f(x)</i>
     */
    public abstract double evaluate(double x);
}

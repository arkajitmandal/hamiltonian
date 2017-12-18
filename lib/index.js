'use strict';
const math = require('mathjs');

// Class for creating Basis Function
class BasisFunction {
  constructor(expression) {
    this.fx = expression;
  }
  // Evaluate the function at a given x
  at(x) {
    return math.eval(this.fx, { x: x });
  }
  // Derivative arround x
  derivative(x, analytical = false) {
    // Numerical Derivative
    if (!analytical) {
      let h = 0.00000001;
      return (this.at(x + h) - this.at(x + h)) / (2 * h);
    }
    // Analytical Derivative
    let Dx = math.eval(math.derivative(this.fx, 'x').toString(), { x: x });
    return Dx;
  }
}

module.exports = {
  BasisFunction: BasisFunction
};

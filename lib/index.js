'use strict';
const math = require('mathjs');
// Const numeric = require('../numeric.js');

// Class for creating Basis Function
class BasisFunction {
  constructor(expression) {
    this.fx = expression;
  }
  // Evaluate the function at a given x
  at(x, fx = this.fx) {
    return math.eval(fx, { x: x });
  }
  // Derivative arround x
  derivative(x, analytical = false) {
    // Numerical Derivative
    if (!analytical) {
      let h = 0.0000001;
      return (this.at(x + h) - this.at(x - h)) / (2 * h);
    }
    // Analytical Derivative
    let Dx = math.eval(math.derivative(this.fx, 'x').toString(), { x: x });
    return Dx;
  }
  // Integrate the function
  // (psi*)(psi)dx
  integrate(from, to, points = 1000) {
    let dx = (to - from) / points;
    let total = 0.0;
    for (var i = 0; i < points; i++) {
      let x = from + i * dx;
      let fx2 = '(' + this.fx + ')^2';
      total += dx / 2 * (this.at(x, fx2) + this.at(x + dx, fx2));
    }
    return total;
  }
  // Normalize the function
  // also return the normalization constant
  normalize(from, to) {
    let N = math.sqrt(this.integrate(from, to));
    this.fx = (1 / N).toString() + '* (' + this.fx + ')';
    return N;
  }
}
module.exports = {
  BasisFunction: BasisFunction
};

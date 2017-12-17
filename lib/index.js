'use strict';
const math = require('mathjs');

// Class for creating Basis Function
class BasisFunction {
  constructor(expression) {
    this.fx = expression;
  }
  // Evaluate the function at a given x
  // eslint-disable-next-line no-unused-vars
  at(x) {
    // eslint-disable-next-line no-eval
    return eval(this.fx);
  }
  // Derivative arround x
  derivative(x, analytical = false) {
    if (!analytical) {
      math.eval(this.fx, { x: x });
    }
  }
}

module.exports = { BasisFunction: BasisFunction };

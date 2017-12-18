'use strict';
const math = require('mathjs');
// Const numeric = require('../numeric.js');

// Class for creating Basis Function
class BasisFunction {
  constructor(expression, limit = [-1, 1]) {
    this.fx = expression;
    this.limit = [Math.min(limit), Math.max(limit)];
    this.upperLimit = Math.max(limit);
    this.lowerLimit = Math.min(limit);
  }
  // Set Limit
  setLimits(limit = [-1, 1]) {
    this.limit = [Math.min(limit), Math.max(limit)];
    this.upperLimit = Math.max(limit);
    this.lowerLimit = Math.min(limit);
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
  // < psi | psi >
  integrate(from = this.lowerLimit, to = this.upperLimit, points = 1000) {
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
  //         ______________
  // 1/N =  V < psi | psi >
  normalize(from = this.lowerLimit, to = this.upperLimit) {
    let N = math.sqrt(this.integrate(from, to));
    this.fx = (1 / N).toString() + '* (' + this.fx + ')';
    return N;
  }
  // Double Derivative
  // d^2/dx^2 (psi)
  doubleDerivative(x, analytical = false) {
    // Numerical 5 point Derivative
    if (!analytical) {
      let h = 0.000001;
      let D1 = this.derivative(x + 2 * h);
      let D2 = this.derivative(x + h);
      let D3 = this.derivative(x - h);
      let D4 = this.derivative(x - 2 * h);
      return (D4 - D1 + 8 * (D2 - D3)) / (12 * h);
    }
    let Dx = math.derivative(this.fx, 'x').toString();
    let D2x = math.derivative(Dx, 'x').toString();
    return math.eval(D2x, { x: x });
  }
  // Kinetic Energy
  // - < psi | d^2/dx^2 | psi >/m
  // eslint-disable-next-line max-params
  kinetic(
    from = this.lowerLimit,
    to = this.upperLimit,
    m = 1.0,
    analytical = false,
    points = 1000
  ) {
    // Value of Kinetic Energy
    let dx = (to - from) / points;
    let total = 0.0;
    for (var i = 0; i < points; i++) {
      let x = from + i * dx;
      let D2x0 = this.doubleDerivative(x, analytical);
      let D2x1 = this.doubleDerivative(x + dx, analytical);
      let fx0 = this.at(x);
      let fx1 = this.at(x + dx);
      total += dx / 2 * (D2x0 * fx0 + D2x1 * fx1);
    }
    return -total / m;
  }
}

module.exports = {
  BasisFunction: BasisFunction
};

const assert = require('assert');
const hamiltonian = require('../index.js');
const math = require('mathjs');

describe('Basis Function', () => {
  it('Creation of a Basis Function', () => {
    let fx = new hamiltonian.BasisFunction('x + 2');
    assert(fx.at(9.0) === 11.0, 'function not created properly');
  });
  it('Derivative of a Basis Function', () => {
    let fx = new hamiltonian.BasisFunction('x^4');
    assert(fx.derivative(0.0, false) === 0.0, 'derivative not 0.0');
    assert(
      fx.derivative(0.0, false) === fx.derivative(0.0, true),
      'analytical and numerical should match at this point'
    );
    assert(
      Math.abs(fx.derivative(5.0, false) - fx.derivative(5.0, true)) < 0.00001,
      ' two derivatives should be close'
    );
  });
  it('Integration of a Basis Function', () => {
    let fx;
    fx = new hamiltonian.BasisFunction('x^2');
    let intfx = math.eval('x^5/5', { x: 2.0 }) - math.eval('x^5/5', { x: -2.0 });
    assert(
      Math.abs(fx.integrate(-2.0, 2.0, 1000) / intfx) > 0.9999,
      'Integration is far from analytical value'
    );
    assert(
      Math.abs(fx.integrate(-2.0, 2.0, 1000) / intfx) < 1.0001,
      'Integration is far from analytical value'
    );
  });
  it('Normalization of a Basis Function', () => {
    let fx;
    fx = new hamiltonian.BasisFunction('x^2 + x^7');
    fx.normalize(-5.0, 5.0);
    assert(
      Math.abs(fx.integrate(-5.0, 5.0, 1000) - 1.0) < 0.00001,
      'Function not normalized'
    );
  });
});

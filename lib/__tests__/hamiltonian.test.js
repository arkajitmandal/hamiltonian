const assert = require('assert');
const hamiltonian = require('../index.js');

describe('Basis Function', () => {
  it('creation of a Basis Function check', () => {
    let fx = new hamiltonian.BasisFunction('x + 2');
    assert(fx.at(9.0) === 11.0, 'function not created properly');
  });
  it('Derivative of a Basis Function check', () => {
    let fx = new hamiltonian.BasisFunction('x^4');
    assert(fx.derivative(0.0, false) === 0.0, 'derivative not 0.0');
    assert(
      fx.derivative(0.0, false) === fx.derivative(0.0, true),
      'analytical and numerical should match at this point'
    );
  });
});

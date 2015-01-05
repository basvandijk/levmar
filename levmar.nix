{ cabal, bindingsLevmar, hmatrix, vector }:

cabal.mkDerivation (self: {
  pname = "levmar";
  version = "HEAD";
  src = ./.;
  buildDepends = [ bindingsLevmar hmatrix vector ];
  meta = {
    homepage = "https://github.com/basvandijk/levmar";
    description = "An implementation of the Levenberg-Marquardt algorithm";
    license = self.stdenv.lib.licenses.bsd3;
    platforms = self.ghc.meta.platforms;
  };
})

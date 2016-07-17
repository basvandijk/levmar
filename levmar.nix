{ mkDerivation, base, bindings-levmar, hmatrix, stdenv, vector }:
mkDerivation {
  pname = "levmar";
  version = "HEAD";
  src = ./.;
  libraryHaskellDepends = [ base bindings-levmar hmatrix vector ];
  homepage = "https://github.com/basvandijk/levmar";
  description = "An implementation of the Levenberg-Marquardt algorithm";
  license = stdenv.lib.licenses.bsd3;
}

let pkgs = import <nixpkgs> {};
in pkgs.haskellPackages.callPackage ./levmar.nix {
     bindingsLevmar = import ../bindings-levmar;
   }

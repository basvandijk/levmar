{ nixpkgs ? import <nixpkgs> {}, compiler ? "default" }:

let

  inherit (nixpkgs) pkgs;

  haskellPackages = if compiler == "default"
                       then pkgs.haskellPackages
                       else pkgs.haskell.packages.${compiler};

  drv = haskellPackages.callPackage (import ./levmar.nix) {
    bindings-levmar = haskellPackages.callPackage (import ../bindings-levmar/bindings-levmar.nix) {};
  };

in

  if pkgs.lib.inNixShell then drv.env else drv

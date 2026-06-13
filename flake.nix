{
  description = "FiniteVarietiesDB TUI";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-26.05";
    poetry2nix = {
      url = "github:nix-community/poetry2nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = { self, nixpkgs, poetry2nix }:
    let
      supportedSystems = [ "x86_64-linux" "aarch64-linux" "x86_64-darwin" "aarch64-darwin" ];
      forAllSystems = nixpkgs.lib.genAttrs supportedSystems;
    in {
      packages = forAllSystems (system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
          p2n = poetry2nix.lib.mkPoetry2Nix { inherit pkgs; };
          
          sharedOverrides = p2n.defaultPoetryOverrides.extend (final: prev: {
            build = pkgs.python312Packages.build;
            pyproject-hooks = pkgs.python312Packages.pyproject-hooks;
          });
        in {
          default = p2n.mkPoetryApplication {
            projectDir = ./.;
            python = pkgs.python312;
            overrides = sharedOverrides;
            meta = {
              license = pkgs.lib.licenses.mit;
            }; 
          };
        });

      devShells = forAllSystems (system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
          p2n = poetry2nix.lib.mkPoetry2Nix { inherit pkgs; };
          
          sharedOverrides = p2n.defaultPoetryOverrides.extend (final: prev: {
            build = pkgs.python312Packages.build;
            pyproject-hooks = pkgs.python312Packages.pyproject-hooks;
          });
        in {
          default = pkgs.mkShell {
            packages = [
              (p2n.mkPoetryEnv { 
                projectDir = ./.; 
                python = pkgs.python312;
                overrides = sharedOverrides;
              })
              pkgs.poetry
            ];
          };
        });
    };
}

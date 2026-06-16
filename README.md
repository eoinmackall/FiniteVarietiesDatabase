# FiniteVarietiesDB.jl

This project aims to provide an open-source database for varities over finite fields.
This project makes critical use of [Oscar.jl](https://github.com/oscar-system/Oscar.jl) (the Open Source Computer Algebra Research system) as a 
foundation for mathematical algorithms and infrastructure and of [DuckDB](https://duckdb.org/) as a foundation for its database
functionality.

## How to use this database:
You can paste the following into a terminal.

```bash
# 1. Clone the repository
git clone https://github.com/eoinmackall/FiniteVarietiesDatabase.git
cd FiniteVarietiesDatabase

# 2. Download the database
wget https://github.com/eoinmackall/FiniteVarietiesDatabase/releases/download/v0.1.0/hypersurfaces_v0.1.0.tar.gz

# 3. Extract the database (merges into data/hypersurfaces/)
tar -xzvf hypersurfaces_v0.1.0.tar.gz
```

The database should now be accessible using infrastructure from DuckDB.

## Using the internal TUI:
A TUI is provided for graphical navigation of the database. The TUI allows for making simple conjunctive database queries.

<div align="center">
  <table>
    <tr>
      <th colspan="3" align="center">Demo</th>
    </tr>
    <tr>
      <td align="center">
        <img width="258" height="282" alt="image" src="https://github.com/user-attachments/assets/126e1828-db9f-4909-91ba-1d364cf42b63" />
      </td>
      <td align="center">
        <img width="258" height="282" alt="image" src="https://github.com/user-attachments/assets/e91fb622-3ead-43f5-81cc-f0999d034564" />
      </td>
      <td align="center">
        <img width="258" height="282" alt="image" src="https://github.com/user-attachments/assets/c77a3241-ada5-4462-99cd-9a0cc2ac8b2a" />
      </td>
    </tr>
  </table>
</div>

#### NixOS
For users on a NixOS system using flakes, you can simply navigate to the root directory and enter:
```bash
nix run .
```

#### Other Linux systems
For users on another Linux distribution, you can install dependencies and run the TUI using [Poetry](https://python-poetry.org/):
```bash
# 1. Install dependencies
poetry install

# 2. launch the TUI locally
poetry run python tui/fvdb.py
```

## TODO:

The following is a list of future plans for this project.

1. Write functionality for transfering the polynomial elements of the database into Oscar(+) objects for use by the user.
2. Add computations of invariants for elements of the database (e.g. zeta function, Picard rank, etc.).
3. Rewrite the hypersurface collection algorithms in C using FLINT.
4. Add functionality for computing lists of isomorphism classes of curves of genus g.

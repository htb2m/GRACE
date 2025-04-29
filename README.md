# GRACE

**GRACE** (Genetic Algorithm for Collagen Engineering) designs heterotrimeric collagen mimetic peptides based on user-defined targets for melting temperature (Tm) and specificity (ΔTm). It can also generate heterotrimers containing user-specified collagen-like sequences.

---

## 1. System Requirements

- **Operating System:**
  - macOS/Linux: Terminal
  - Windows: Command Prompt, PowerShell, or Windows Terminal
- **Compiler:**
  - C++11 or later

---

## 2. Compilation
- Download the source code
- Navigate to the directory containing the source files: `main.cpp`  `Functions.cpp`  `SCEPTTr.cpp`  `Functions.hpp`  `SCEPTTr.hpp`


	```
 	cd /path/to/dir
	```

- Compile using `g++`:

	```
 	g++ -std=c++11 -o GRACE main.cpp Functions.cpp SCEPTTr.cpp
	```

This will generate an executable named `GRACE`.

---

## 3. Required Files for Execution

- **`AminoAcids.csv`**  
  Each cell of the table is either `1` or `0`. Set to `0` to exclude specific amino acids at Xaa or Yaa positions.

- **`parameters.txt`**  
  Configuration file containing scoring parameters.

Ensure these files are in the same directory as the compiled executable.

---

## 4. Running GRACE

- To execute:

	```
 	./GRACE
 	```

- Follow the on-screen prompts to input search conditions.

## 5. Notes
- Search time may range from minutes to hours depending on sequence length and conditions.

- Please expect that there are chances that the algorithm cannot find a set of sequences that satisfy all given conditions. In that case, the algorithm is terminated after 500,000 iterations and outputs the closest results found. 

- Please expect instance when the algorithm shows minimal progress over several generations. In that case, re-initializing the population by starting the algorithm over can help escape the trap by setting the search at a different location.

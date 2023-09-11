
# Global Similarity and Distance Computer

This project computes the affine indel gap model global sequence
alignment, using the provided sequence strings and score values.






## Authors

- [@carnol29](https://www.github.com/carnol29) / [@ArnoCame](https://www.github.com/ArnoCame)
- [@coom-slayer](https://www.github.com/coom-slayer)


## Running Tests

To run a test, run the following command:
## Run Locally

Clone the project

```bash
  git clone https://github.com/carnol29/global-alignment
```

Go to the project directory

```bash
  cd global-alignment
```

Run the program

```bash
  python3 main.py ["similarity"|"distance"] {blosum_file} {sequences_file}
```


## Usage/Examples

Running a global distance alignment with test files:

```bash
python3 main.py distance input_files/blosum62.txt input_files/sequences.txt
```

Running a global similarity alignment with test files:

```bash
python3 main.py similarity input_files/blosum62.txt input_files/sequences.txt
```

### Instruction to replicate the results presented in the paper

0. Make sure the is an internet connection and using Linux/Unix/macOS/WSL with amd64 architecture (for gobnilp)
1. Download Docker from docker.com
2. Head to the same directory as this file and type 
    ```sh
        docker build -t dnc .  
    ```
3. Start a docker container by
    ```sh
        docker run -it --privileged -w /neurips24 -v $(pwd):/neurips24 dnc bash
    ```
4. To run the 1000 random seeds (starting at 1) simulations using 50 cores in parallel (20 for each core) do:
    ```sh
        ./runsims.sh 1 50 20 results/ tmpfile.csv
    ```
5. Join the results into one CSV file called paper_joint.csv by
    ```sh
        Rscript R/benchmark_opruner.R --seeds_from 1 --seeds_to 1000 --output_dir results/ --filename paper_joint.csv
    ```
6. To generate all figures except for Figure 4 run 
    ```sh
        Rscript R/paper_benchmarks.R
    ```
7. To generate Figure 4 do
    ```sh
        Rscript R/gobnilp_comparison.R
    ```

The plots are saved in a subfolder called figures/
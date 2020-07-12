#Qualitative Switch (qu_switch)

Wright-Fisher simulations with mutation for two-locus three-allele system.

##Syntax
```python
python simulate.py <POPULATION_SIZE> <TOTAL_GENERATIONS> <TOTAL_REPLICATES> <BLOCK_ID>
```

##Example
To simulate 5 replicates of 1000 individuals for 300 generations.
```python
python simulate.py 1000 300 5 1
```

##Parallelization
On *Windows* open as many command prompts as number of parallel jobs and run with different block IDs:
Command prompt 1
```
python simulate.py 1000 300 5 1
```
Command prompt 2
```
python simulate.py 1000 300 5 2
```

On *Linux* or *macOS* using `parallel`
```shell
parallel -k --lb python simulate.py ::: 1000 ::: 300 ::: 5 ::: 1 2
```
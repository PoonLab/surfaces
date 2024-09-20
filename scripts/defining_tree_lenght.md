# Determining number of sequences and time scale

## 1. Simulate a random tree with 100 tips:
```R
require ape
r.tree <- rtree(100, rooted = TRUE,
                tip.label = NULL)
write.tree(r.tree, file = "random_tree_100.nwk", append = FALSE,
           digits = 10, tree.names = FALSE)

```
## 2. Simulate sequences with INDELible
Use `create_control_indelible.py` to create a `control.txt` file:
```
python3 create_control_indelible.py random_tree_100.nwk -o random_sim -r -s 0.5
```
Note: when using the option `-r`, the script runs indelible from the control file created. 
Run output will be stored at `random_sim`

Alternatively, run indelible from terminal:
```
indelible control.txt
```
## 3. Measure selection with HyPhy
```
hyphy
1
8
output_file.fas  # sequences from indelible
random_tree100.nwk  # tree from ape
```

## Alternativelyu: run 2 and 3 in a single script for multiple tree lengths
To run indelible and measure selection, all at once, use the script `indelible_plus_selection.py`.
This script uses the function `run_selection_pipeline` from `selection_by_cluster.py`.
```
python3 indelible_plus_selection.py ../data_find_tree_length/random_tree_100.nwk -p ../data_find_tree_length/ -r 1
```

## 4. Compare selection measurments with simulation inputs
Simulation inputs for hyphy were:
```
omegas = [
0.0510204, 0.1530612, 0.2551020, 0.3571428, 0.4591836, 0.5612244, 0.6632652, 0.7653060, 0.8673468,
0.9693876, 1.0714284, 1.1734692, 1.2755100, 1.3775508, 1.4795916, 1.5816324, 1.6836732, 1.7857140,
1.8877548, 1.9897956, 2.0918364, 2.1938772, 2.2959180, 2.3979588, 2.4999996, 2.6020404, 2.7040812,
2.8061220, 2.9081628, 3.0102036, 3.1122444, 3.2142852, 3.3163260, 3.4183668, 3.5204076, 3.6224484,
3.7244892, 3.8265300, 3.9285708, 4.0306116, 4.1326524, 4.2346932, 4.3367340, 4.4387748, 4.5408156,
4.6428564, 4.7448972, 4.8469380, 4.9489788, 5.0510196]

prop = [
    1.063764e-01, 1.464867e-01, 1.401630e-01, 1.223917e-01, 1.023007e-01, 8.333079e-02, 6.673162e-02, 
    5.279576e-02, 4.139382e-02, 3.222714e-02, 2.495008e-02, 1.922789e-02, 1.476163e-02, 1.129629e-02,
    8.620596e-03, 6.562952e-03, 4.985996e-03, 3.780965e-03, 2.862470e-03, 2.163923e-03, 1.633688e-03,
    1.231906e-03, 9.279287e-04, 6.982665e-04, 5.249686e-04, 3.943506e-04, 2.960034e-04, 2.220245e-04,
    1.664245e-04, 1.246710e-04, 9.333909e-05, 6.984377e-05, 5.223631e-05, 3.904912e-05, 2.917804e-05,
    2.179306e-05, 1.627075e-05, 1.214322e-05, 9.059520e-06, 6.756626e-06, 5.037501e-06, 3.754636e-06,
    2.797655e-06, 2.084012e-06, 1.551998e-06, 1.155507e-06, 8.600995e-07, 6.400653e-07, 4.762155e-07]
```

If `[printrates] TRUE`, it will print a `RATES.txt` file like:
```
Site	Class	Partition	Inserted?
1	3	1
2	3	1
3	9	1
4	8	1
5	4	1
6	4	1
7	16	1
8	6	1
9	2	1
10	4	1
11	0	1
12	9	1
```

Now I need to map each site to the omega value it corresponds to, and then I can compare site by site against fubar output

**Note:** HyPhy doesn't like INDELible output when `indel rate = 0`, so I re-aligned the file with `mafft` and now it runs properly.

**Note:** To import the `RATES.txt` files to R from the INDELible outpit, we need to remove the first 9 lines with information about the simulation:
```
sed -i '1,9d' *_RATES.txt 
``` 

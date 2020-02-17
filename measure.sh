make clean
make `OPTIONS.TXT`

for n in 256 257 329 511 512 513 631 768 843 960 1023 1024 1025 1374 1536 1720 1845 2047 2048 2049
do
r=$[10+78125/($n/120)/($n/120)/($n/120)]
echo "n=$n, r=$r"
./mmpy -n $n -r $r | tee results/data_$n.txt | grep 'Device computation time'
done

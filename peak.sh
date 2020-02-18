mkdir -p peak

for n in 256 512 1024 2048
do
for bx in 4 8 16 32
do
for by in 4 8 16 32
do
r=$[(10+200000/($n/120+1)/($n/120+1)/($n/120+1))*$bx*$by/16/16]
if [[ $[$bx*$by] -le $[16*16] ]]
then
make clean >/dev/null 2>&1
make "bx=$bx" "by=$by" >/dev/null 2>&1
echo -n "n=$n, r=$r, bx=$bx, by=$by, " | tee -a "peak/summary.txt"
./mmpy -n $n -r $r -R | tee "peak/result_${n}_${bx}_${by}.txt" | grep 'Device computation time' | tee -a "peak/summary.txt"
fi
done
done
done

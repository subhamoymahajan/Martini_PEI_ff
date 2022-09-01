rm -f performance.dat
time=0
for i in {1..8}
do
	a=`grep 'Time:' CG${i}/npt.log | awk '{print$3}'`
	time=$(bc -l <<<"${time}+${a}")
	echo -e "CG${i}\t NPT:\t $a s" >> performance.dat
	a=`grep 'Time:' CG${i}/md_1.log | awk '{print$3}'`
	time=$(bc -l <<<"${time}+${a}")
	echo -e "CG${i}\t MD:\t $a s" >> performance.dat
	if [ $i -lt 8 ]
	then
		a=`grep 'Time:' CG${i}_th/npt.log | awk '{print$3}'`
		time=$(bc -l <<<"${time}+${a}")
		echo "CG${i}_th\t NPT:\t $a s" >> performance.dat
		a=`grep 'Time:' CG${i}_th/md_1.log | awk '{print$3}'`
		time=$(bc -l <<<"${time}+${a}")
		echo "CG${i}_th\t MD:\t $a s" >> performance.dat
		a=`grep 'Time:' CG${i}_K/npt.log | awk '{print$3}'`
		time=$(bc -l <<<"${time}+${a}")
		echo "CG${i}_K\t NPT:\t $a s" >> performance.dat
		a=`grep 'Time:' CG${i}_K/md_1.log | awk '{print$3}'`
		time=$(bc -l <<<"${time}+${a}")
		echo "CG${i}_K\t MD:\t $a s" >> performance.dat
	fi

done
time=$(bc -l <<<"${time}/60.0/60.0")
echo -e "Total time: $(printf %0.2f $time) hours\n$(cat performance.dat)" > performance1.dat
rm performance.dat
mv performance1.dat performance.dat

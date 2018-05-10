rm now.csv
echo 'score' >> now.csv
#cd ../vis
for i in {1..1000}
do
    run1=$(java -jar tester.jar -exec "./a" -seed $i -novis)
    echo $run1 >> ../tester/now.csv
done

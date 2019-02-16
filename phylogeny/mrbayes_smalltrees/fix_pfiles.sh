#This should fix all parameter files in a folder, deleting everything after incomplete lines:

for test_file in *.p
do

	#get number o tabs a line should have
	tab_num=$(awk -F'|' 'BEGIN{print "count", "lineNum"}{print gsub(/\t/,"") "\t" NR}' $test_file | head -n 3 | tail -n 1 | cut -f 1)

	#find if any line does not have this number of tabs - a sign of problems
	problem_line=$(awk -F'|' 'BEGIN{print "count", "lineNum"}{print gsub(/\t/,"") "\t" NR}' $test_file | tail -n +3 | grep -v -E "^$tab_num\s" | head -n 1 | cut -f 2)

	#if problem line, remove everything from that line on:
	head -n $(($problem_line-1)) $test_file > tmp && mv tmp $test_file
        rm tmp
done

for test_file in *.t
do

        #get number of characters a line should have
        par_num=$(awk -F'|' 'BEGIN{print "count", "lineNum"}{print gsub(/)/,"") "\t" NR}' $test_file | grep -v -E "^0\s" | head -n 3 | tail -n 1 | cut -f 1) 

        #find if any line does not have this number of tabs - a sign of problems
        problem_line=$(awk -F'|' 'BEGIN{print "count", "lineNum"}{print gsub(/)/,"") "\t" NR}' $test_file | tail -n +2 | grep -v -E "^0\s" | grep -v -E "^$par_num\s" | head -n 1 | cut -f 2)  
    
        #if problem line, remove everything from that line on:
        head -n $(($problem_line-1)) $test_file > tmp && mv tmp $test_file
        rm tmp
done


for test_file in *.mcmc
do

	#get number o tabs a line should have
	tab_num=$(awk -F'|' 'BEGIN{print "count", "lineNum"}{print gsub(/\t/,"") "\t" NR}' $test_file | head -n 10 | tail -n 1 | cut -f 1)

	#find if any line does not have this number of tabs - a sign of problems
	problem_line=$(awk -F'|' 'BEGIN{print "count", "lineNum"}{print gsub(/\t/,"") "\t" NR}' $test_file | tail -n +10 | grep -v -E "^$tab_num\s" | head -n 1 | cut -f 2)

	#if problem line, remove everything from that line on:
	head -n $(($problem_line-1)) $test_file > tmp && mv tmp $test_file
        rm tmp
done


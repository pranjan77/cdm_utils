#conda env export > environment.yml
sed 's/\(.*\)=.*/\1/' $1

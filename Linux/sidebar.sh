for i in `ls -tr | tr "\t" "\n"`;
do
i=${i%.md*};
echo - [${i}]\(/题目/${i}.md\) >>../_sidebar.md
done
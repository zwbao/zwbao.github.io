for i in `ls -tr | tr "\t" "\n"`;
do
i=${i%.md*};
echo - [${i}]\(/Python/基础/${i}.md\) >>../_sidebar.md
done
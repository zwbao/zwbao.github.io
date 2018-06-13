for i in `ls -tr | tr "\t" "\n"`;
do
i=${i%.md*};
echo - [${i}]\(/R/${i}.md\) >>_sidebar.md
done
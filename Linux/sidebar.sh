for i in *.md;
do
i=${i%.md*};
echo - [${i}]\(/Linux基础/${i}.md\) >>_sidebar.md
done
if [ "$#" -ne 1 ]
then
   echo "You should pass as argument the txt file name"
   exit
fi
if [ ! -e $1 ]
then
   echo "File" $1 "does not exist"
   exit
fi

echo processing file $1
ls -l $1

root -l -q ConvertTxtToRoot.C+\(\"$1\"\)

rootFileName=$1
rootFileName=${rootFileName/.txt/.root}
ls -l $rootFileName

root -l -q ProcessFile.C+\(\"$rootFileName\"\)

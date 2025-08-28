export CVSROOT=/afs/cern.ch/ams/Offline/CVS
cd /afs/cern.ch/ams/Offline/CVS/AMS/CC
for file in *.C,v *.h,v *.cc,v *.cpp,v; do
    if [ -f "$file" ]; then
        filename="${file%,v}"
        cvs co -p "AMS/CC/$filename" 2>/dev/null | grep -A 2 -B 2 -n "getTrackEmissionPoint" | sed "s/^/$filename:/"
    fi
done > /afs/cern.ch/work/z/zuhao/public/yanzx/IsoAnalys/richpos.log
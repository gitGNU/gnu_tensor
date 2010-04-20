#!/bin/sh

for i in char float int long short; do
    cp tensor_TYPE.h tensor_${i}.h
    perl -p -i.bak -e "s/NAME/${i}/g" tensor_${i}.h
    perl -p -i.bak -e "s/TYPE/${i}/g" tensor_${i}.h
done



cp tensor_TYPE.h tensor_long_double.h
perl -p -i.bak -e "s/NAME/long_double/g" tensor_long_double.h
perl -p -i.bak -e "s/TYPE/long double/g" tensor_long_double.h

cp tensor_TYPE.h tensor_uchar.h
perl -p -i.bak -e "s/NAME/uchar/g" tensor_uchar.h
perl -p -i.bak -e "s/TYPE/unsigned char/g" tensor_uchar.h

cp tensor_TYPE.h tensor_uint.h
perl -p -i.bak -e "s/NAME/uint/g" tensor_uint.h
perl -p -i.bak -e "s/TYPE/unsigned int/g" tensor_uint.h

cp tensor_TYPE.h tensor_ulong.h
perl -p -i.bak -e "s/NAME/ulong/g" tensor_ulong.h
perl -p -i.bak -e "s/TYPE/unsigned long/g" tensor_ulong.h

cp tensor_TYPE.h tensor_ushort.h
perl -p -i.bak -e "s/NAME/ushort/g" tensor_ushort.h
perl -p -i.bak -e "s/TYPE/unsigned short/g" tensor_ushort.h


rm *.bak

# To delete:
# rm `\ls -1 tensor_*.h | grep -v double.h | grep -v complex | grep -v TYPE | grep -v utili`

BIN_DIR=$(pwd)/bin
LIB_DIR=$(pwd)/lib
SRC_DIR=$(pwd)/src
INCLUDE_DIR=$(pwd)/include
if [ ! -d "$BIN_DIR" ]; then
  	mkdir "$BIN_DIR"
fi
/usr/bin/g++ -L$LIB_DIR -Wl,-rpath=$LIB_DIR --std=c++11 -g -O2 $SRC_DIR/PanDepth.cpp -lhts -lpthread -lz -ldeflate -I $INCLUDE_DIR -o $BIN_DIR/pandepth
if [ $? -ne 0 ]; then
  echo "Installation of PanDepth has failed."
else
  echo "The installation of PanDepth has been successfully completed."
fi

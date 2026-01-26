CC = g++
CFLAGS = -Wall -I/usr/include/ -I/usr/include/eigen3 -I/usr/include/eigen3 -L/usr/local/lib
LIBS = -lgsl -lgslcblas -lm -larmadillo

BIN_DIR = bin
SRC_DIR = src

# Executable targets
TARGETS = \
  $(BIN_DIR)/anl \
  $(BIN_DIR)/sgl \
  $(BIN_DIR)/num \
  $(BIN_DIR)/fro 

all: $(TARGETS)

$(BIN_DIR)/anl: $(SRC_DIR)/anl_time.cxx
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

$(BIN_DIR)/sgl: $(SRC_DIR)/single.cxx
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

$(BIN_DIR)/num: $(SRC_DIR)/num_time.cxx
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

$(BIN_DIR)/fro: $(SRC_DIR)/frohlich.cxx
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

clean:
	rm -f $(BIN_DIR)/* $(SRC_DIR)/*.o

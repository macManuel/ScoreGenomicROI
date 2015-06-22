TARGET = ScoreGenomicROI
PREFIX =  /usr/local

HEADER_PATH = $(PREFIX)/include/
LIB_PATH = $(PREFIX)/lib/

HEADER = ScoreGenomicROI.h
OBJECTS = ScoreGenomicROI.o
LIBS = libFastScience.dylib

CC = g++
CFLAGS = -std=c++11 -Wall -I$(HEADER_PATH)

$(TARGET): $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) -o $(TARGET) $(LIB_PATH)/$(LIBS) $(OBJECTS) main.cpp

ScoreGenomicROI.o : ScoreGenomicROI.cpp ScoreGenomicROI.h
	$(CC) $(CFLAGS) -c -o ScoreGenomicROI.o ScoreGenomicROI.cpp 

install:
	cp -p $(TARGET) $(PREFIX)/bin/

clean:
	rm -r $(OBJECTS)
	rm -r $(TARGET)
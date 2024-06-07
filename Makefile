LIBS_ROOT = ../../pacs-examples/Examples

CXX = mpicxx
CXXFLAGS = -fopenmp -O3 -std=c++20
CPPFLAGS = -I${LIBS_ROOT}/include
LDFLAGS = -L${LIBS_ROOT}/lib
LIBS = -lmuparser

ifdef _SCALABILITY
    CXXFLAGS += -D_SCALABILITY
endif

TARGET = laplace_solver
INCLUDE_DIR = include
SRC_DIR = src
SRCS = $(SRC_DIR)/main.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I$(INCLUDE_DIR) -o $@ $^ $(LDFLAGS) $(LIBS)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I$(INCLUDE_DIR) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
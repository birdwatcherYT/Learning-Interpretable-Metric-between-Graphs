###################################################
# Makefile (mylib.a 利用する場合)
# コンパイルしたいソースがあるフォルダと同じところに置く
# ディレクトリ階層図
# this
# |--- *.cpp          (コンパイルしたいソースファイル)
# |--- Makefile
# |--- run.exe        (生成される実行ファイル)
# |--- obj            (生成されるフォルダ)
# |    |--- *.o       (生成されるオブジェクトファイル)
# |    |--- *.d       (生成される依存関係ファイル)
# |
# |--- lib            (用意するフォルダ)
#      |--- *cpp      (ライブラリのソースファイル)
#      |--- *hpp      (ライブラリのヘッダーファイル)
#      |--- mylib.a   (生成されるライブラリ)
#      |--- obj       (生成されるフォルダ)
#           |--- *.o  (生成されるオブジェクトファイル)
#           |--- *.d  (生成される依存関係ファイル)
#
# 実行コマンド
# make         (更新されたファイルのみコンパイル)
# make all     (クリーンしてビルド)
# make clean   (クリーン)
###################################################
AR       = ar
ARFLAGS  = rus

CXX      = g++
CXXFLAGS = -g -MMD -MP -Wall -Wextra -Winit-self -Wno-unused-parameter -O3 -fopenmp

RM       = rm -f
LDFLAGS  = -fopenmp
INCLUDE  = 
#実行ファイル
TARGET   = ./run.exe
#オブジェクトファイル保存フォルダ
OBJDIR   = ./obj
#ソースファイル
SOURCES  = $(wildcard *.cpp)
#オブジェクトファイル
OBJECTS  = $(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))
#依存関係ファイル
DEPENDS  = $(OBJECTS:.o=.d)

#ライブラリのあるフォルダ
LIB_DIR     = 
#ライブラリの相対パス
LIBS        = $(LIB_DIR)
#ライブラリのヘッダーファイル
LIB_HPP     = $(wildcard $(LIB_DIR)/*.hpp)
#ライブラリのオブジェクトファイル保存フォルダ
LIB_OBJ_DIR = $(LIB_DIR)/obj
#ライブラリのソースファイル
LIB_SRC     = $(wildcard $(LIB_DIR)/*.cpp)
#ライブラリのオブジェクトファイル
LIB_OBJ     = $(subst $(LIB_DIR),$(LIB_OBJ_DIR),$(LIB_SRC:.cpp=.o))
#ライブラリの依存関係ファイル
LIB_DEP     = $(LIB_OBJ:.o=.d)

#run.exeの作成
$(TARGET): $(OBJECTS) $(LIBS)
	$(CXX) -o $@ $^ $(LDFLAGS)

#ライブラリ作成
$(LIBS): $(LIB_OBJ)
	$(AR) $(ARFLAGS) $@ $^

#libフォルダ内にobjフォルダ作成後、libフォルダ内の*.cppをコンパイルし、*.oの生成
$(LIB_OBJ_DIR)/%.o: $(LIB_DIR)/%.cpp
	@if [ ! -d $(LIB_OBJ_DIR) ];\
	then echo "mkdir -p $(LIB_OBJ_DIR)";mkdir -p $(LIB_OBJ_DIR);\
	fi
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

#コンパイルしたいソースファイルのコンパイル
$(OBJDIR)/%.o: %.cpp $(LIB_HPP)
	@if [ ! -d $(OBJDIR) ];\
	then echo "mkdir -p $(OBJDIR)";mkdir -p $(OBJDIR);\
	fi
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

#クリーンしてビルド
.PHONY: all
all: clean $(TARGET)

#クリーン
.PHONY:clean
clean:
	$(RM) $(OBJECTS) $(DEPENDS) $(TARGET) $(LIB_OBJ) $(LIB_DEP)

-include $(DEPENDS)

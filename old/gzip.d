import std.zlib;
import std.stdio;
import std.range;
import std.traits;

class GzipInputRange {
  UnCompress uncompressObj;
  File f;
  auto CHUNKSIZE = 0x4000;
  ReturnType!(f.byChunk) chunkRange;
  bool exhausted;
  char[] uncompressedBuffer;
  size_t bufferIndex;
  
  this(string filename) {
    f = File(filename, "r");
    chunkRange = f.byChunk(CHUNKSIZE);
    uncompressObj = new UnCompress(HeaderFormat.gzip);
    load();
  }
  
  void load() {
    if(!chunkRange.empty) {
      auto raw = chunkRange.front.dup;
      chunkRange.popFront();
      uncompressedBuffer = cast(char[])uncompressObj.uncompress(raw);
      // uncompressedBuffer = cast(char[])(uncompressObj.uncompress(raw).dup);
      
      bufferIndex = 0;
    }
    else {
      if(!exhausted) {
        uncompressedBuffer = cast(char[])uncompressObj.flush();
        // uncompressedBuffer = cast(char[])(uncompressObj.flush().dup);
        exhausted = true;
        bufferIndex = 0;
      }
      else
        uncompressedBuffer.length = 0;
    }
  }
  
  @property char front() {
    return uncompressedBuffer[bufferIndex];
  }
  
  void popFront() {
    bufferIndex += 1;
    if(bufferIndex >= uncompressedBuffer.length) {
      load();
      bufferIndex = 0;
    }
  }
  
  @property bool empty() {
    return uncompressedBuffer.length == 0;
  }
}

class GzipByLine {
  GzipInputRange range;
  char[] buf;
  
  this(string filename) {
    this.range = new GzipInputRange(filename);
    popFront();
  }
  
  @property bool empty() {
    return buf.length == 0;
  }
  
  void popFront() {
    buf.length = 0;
    while(!range.empty && range.front != '\n') {
      buf ~= range.front;
      range.popFront();
    }
    range.popFront();
  }
  
  string front() {
    return buf.idup;
  }
}

class GzipOut {
  Compress compressObj;
  File f;
    
  this(string filename) {
    f = File(filename, "w");
    compressObj = new Compress(HeaderFormat.gzip);
  }
  
  void compress(string s) {
    try {
      auto compressed = compressObj.compress(s.dup);
      f.rawWrite(compressed);
    }
    catch (Exception e) {
      stderr.writeln("trying to compress: ", s);
      throw e;
    }
  }

  void finish() {
    auto compressed = compressObj.flush();
    f.rawWrite(compressed);
  }
}


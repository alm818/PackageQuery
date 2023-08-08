#pragma once

class AbstractWalker {
   public:
      virtual int step() = 0;
      virtual ~AbstractWalker(){
      }
};
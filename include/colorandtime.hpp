#ifndef GAIA_COLORS_HPP
#define GAIA_COLORS_HPP

#include <ostream>
#include <sys/time.h>

namespace chronos{
  inline double time(){
    struct timeval tmp_time;
    gettimeofday(&tmp_time,NULL);
    return tmp_time.tv_sec+(tmp_time.tv_usec*1.0e-6L);
  }
} // CHRONOS

namespace colors{
  enum Code {
    FG_RED        = 31,
    FG_GREEN      = 32,
    FG_BLUE       = 34,
    FG_DEFAULT    = 39,
    FG_BLACK      = 30,
    FG_YELLOW     = 33,
    FG_MAGENTA    = 35,
    FG_CYAN       = 36,
    FG_LIGHTGREY  = 37,
    FG_DARKGREY   = 90,
    FG_WHITE      = 97,
    BG_RED        = 41,
    BG_GREEN      = 42,
    BG_BLUE       = 44,
    BG_DEFAULT    = 49,
    TXT_BOLD      = 1,
    TXT_DIM       = 2,
    TXT_UNDERLINE = 4,
    TXT_BLINK     = 5,
    TXT_DEFAULT   = 0
  };
  class metamorphe {
      Code code;
  public:
    metamorphe(Code pCode) : code(pCode) {}
    friend std::ostream&
    operator<<(std::ostream& os, const metamorphe& mod) {
      return os << "\033[" << mod.code << "m";
    }
  };
} // COLORS

colors::metamorphe f_red      (colors::FG_RED       );
colors::metamorphe f_blue     (colors::FG_BLUE      );
colors::metamorphe f_green    (colors::FG_GREEN     );
colors::metamorphe f_def      (colors::FG_DEFAULT   );
colors::metamorphe f_black    (colors::FG_BLACK     );
colors::metamorphe f_yellow   (colors::FG_YELLOW    );
colors::metamorphe f_magenta  (colors::FG_MAGENTA   );
colors::metamorphe f_cyan     (colors::FG_CYAN      );
colors::metamorphe f_lgrey    (colors::FG_LIGHTGREY );
colors::metamorphe f_dgrey    (colors::FG_DARKGREY  );
colors::metamorphe f_white    (colors::FG_WHITE     );
colors::metamorphe b_red      (colors::BG_RED       );
colors::metamorphe b_blue     (colors::BG_BLUE      );
colors::metamorphe b_green    (colors::BG_GREEN     );
colors::metamorphe b_def      (colors::BG_DEFAULT   );
colors::metamorphe t_bold     (colors::TXT_BOLD     );
colors::metamorphe t_dim      (colors::TXT_DIM      );
colors::metamorphe t_underline(colors::TXT_UNDERLINE);
colors::metamorphe t_blink    (colors::TXT_BLINK    );
colors::metamorphe t_def      (colors::TXT_DEFAULT  );

#endif

/*
 *  Copyright (C) 2018   Malte Brunn
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <stdio.h>
#include <map>
#include <string>
#include <functional>
#include <sstream>
#include <tuple>
//--------------------------------------------------------------------------------------------------
#ifndef _ARGVPARSER_HPP
#define _ARGVPARSER_HPP
//--------------------------------------------------------------------------------------------------
class ARGVParser {
private:
  typedef std::function<int(int, char **)> callback_functor;
  std::map<std::string, std::tuple<int, int, callback_functor> > _callbacks;
  std::map<std::string, std::string> _aliases;
public:
  ARGVParser () {
    _callbacks.clear();
    _aliases.clear();
  }
  void alias(std::string param, std::string alias) {
    if (_callbacks.count(alias) == 0)
      _aliases[alias] = param;
  }
  template <typename Type>
  void bind (std::string param, Type &var) {
    _aliases.erase(param);
    _callbacks[param] = std::make_tuple(1, 0, [&var](int ac, char **av) -> int {
      if (ac != 1) return 0;
      std::istringstream stream(av[0]);
      Type tmp;
      stream >> tmp;
      if (stream.fail()) return 0;
      var = tmp;
      return 1;
    });
  }
  template <typename Type, typename Functor>
  void bind (std::string param, Type &var, Functor lambda) {
    _aliases.erase(param);
    _callbacks[param] = std::make_tuple(1, 0, [&var, lambda](int ac, char **av) -> int {
      if (ac != 1) return 0;
      try {
        Type tmp;
        tmp = lambda(av[0]);
        var = tmp;
      } catch (...) {
        return 0;
      }
      return 1;
    });
  }
  void bind (std::string param, std::function<int(int, char **)> lambda, int cnt=1) {
    _callbacks[param] = std::make_tuple(cnt, 0, lambda);
  }
  void bind (std::string param, std::function<void()> lambda) {
    _callbacks[param] = std::make_tuple(0, 0, [lambda](int, char**) -> int {
      lambda();
      return 0;
    });
  }
  void bind (std::string param) {
    _callbacks[param] = std::make_tuple(0, 0, [](int, char**) -> int { return 0; });
  }
  void exec (int argc, char **argv) {
    for (int i = 0; i < argc; ++i) {
      std::string param = argv[i];
      if (_aliases.count(param)) {
        auto alias = _aliases[param];
        param = alias;
      }
      if (_callbacks.count(param)) {
        auto &functor = _callbacks[param];
        int cnt = std::get<0>(functor);
        if (cnt < 0 || i + cnt >= argc) cnt = argc - i - 1;
        int used = std::get<2>(functor)(cnt, &argv[i+1]);
        if (used || cnt == 0) {
          std::get<1>(functor) += 1;
        }
        i += used;
      }
    }
  }
  int count (std::string param) {
    if (_aliases.count(param)) {
      auto alias = _aliases[param];
      param = alias;
    }
    if (_callbacks.count(param)) {
      return std::get<1>(_callbacks[param]);
    }
    return 0;
  }
  void reset () {
    for (auto &callback : _callbacks) {
      std::get<1>(callback.second) = 0;
    }
  }
  void help () {
    for (auto &callback : _callbacks) {
      printf("[%s", callback.first.c_str());
      for (auto &alias : _aliases)
        if (alias.second == callback.first)
          printf(",%s",alias.first.c_str());
      printf("]");
    }
    printf("\n");
  }
};
//--------------------------------------------------------------------------------------------------
#endif // _ARGVPARSER_HPP

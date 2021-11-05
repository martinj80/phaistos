// program_options_short_format_parser.h --- Parser for command line options in short format
// Copyright (C) 2006-2013 Wouter Boomsma
//
// This file is part of Phaistos 
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.

// Inspired by boost::property_tree::json_parser.hpp


#ifndef PROGRAM_OPTIONS_SHORT_FORMAT_PARSER_H
#define PROGRAM_OPTIONS_SHORT_FORMAT_PARSER_H

#include <boost/property_tree/ptree.hpp>
#include <boost/spirit/include/classic.hpp>
#include <boost/property_tree/detail/file_parser_error.hpp>

namespace phaistos {
namespace option_parser {

template<class Ptree>
void write_short_option_format_helper(std::basic_ostream<typename Ptree::key_type::value_type> &stream,
                                      const Ptree &pt, int recursion_level) {

     typedef typename Ptree::key_type::value_type Ch;
     typedef typename std::basic_string<Ch> Str;

     // Value or object or array
     if (recursion_level > 0 && pt.empty()) {
          // Write value
          Str data = pt.template get_value<Str>();
          stream << data;

     }
     else {
          // Write object
          if (recursion_level > 0)
               stream << Ch('[');
          typename Ptree::const_iterator it = pt.begin();
          for (; it != pt.end(); ++it)
          {
               stream << it->first << Ch(':');
               write_short_option_format_helper(stream, it->second, recursion_level + 1);
               if (boost::next(it) != pt.end())
                    stream << Ch(',');
          }
          if (recursion_level > 0)
               stream << Ch(']');
     }
}

template<class Ptree>
void write_short_option_format(std::basic_ostream<typename Ptree::key_type::value_type> &stream,
                               const Ptree &pt) {
     
     write_short_option_format_helper(stream, pt, 0);
     
}




//! Parser context (defines functionality at nodes during parsing)
template<class Ptree>
struct context
{

     typedef typename Ptree::key_type::value_type Ch;
     typedef std::basic_string<Ch> Str;
     typedef typename std::vector<Ch>::iterator It;
        
     Str string;
     Str name;
     Ptree root;
     std::vector<Ptree *> stack;

     struct a_object_s {
          context &c;
          a_object_s(context &c): c(c) { }
          void operator()(It, It) const {
               if (c.stack.empty())
                    c.stack.push_back(&c.root);
               else
               {
                    Ptree *parent = c.stack.back();
                    Ptree *child = &parent->push_back(std::make_pair(c.name, Ptree()))->second;
                    c.stack.push_back(child);
                    c.name.clear();
               }
          }
     };
        
     struct a_object_e {
          context &c;
          a_object_e(context &c): c(c) { }
          void operator()(It, It) const {
               BOOST_ASSERT(c.stack.size() >= 1);
               c.stack.pop_back();
          }
     };

     struct a_name {
          context &c;
          a_name(context &c): c(c) { }
          void operator()(It, It) const {
               c.name.swap(c.string);
               c.string.clear();
          }
     };

     struct a_string_val {
          context &c;
          a_string_val(context &c): c(c) { }
          void operator()(It, It) const {
               BOOST_ASSERT(c.stack.size() >= 1);
               c.stack.back()->push_back(std::make_pair(c.name, Ptree(c.string)));
               c.name.clear();
               c.string.clear();
          }
     };

     struct a_literal_val {
          context &c;
          a_literal_val(context &c): c(c) { }
          void operator()(It b, It e) const {
               BOOST_ASSERT(c.stack.size() >= 1);
               c.stack.back()->push_back(std::make_pair(c.name,
                                                        Ptree(Str(b, e))));
               c.name.clear();
               c.string.clear();
          }
     };

     struct a_char {
          context &c;
          a_char(context &c): c(c) { }
          void operator()(It b, It e) const {
               c.string += *b;
          }
     };

};

// Grammar
template<class Ptree>
struct short_option_grammar :
          public boost::spirit::classic::grammar<short_option_grammar<Ptree> > {
        
     typedef context<Ptree> Context;
     typedef typename Ptree::key_type::value_type Ch;

     mutable Context c;
        
     template<class Scanner>
     struct definition {
          
          boost::spirit::classic::rule<Scanner> root, object, member, value, symbol, string, number;
          boost::spirit::classic::rule<typename boost::spirit::classic::lexeme_scanner<Scanner>::type> character;

          definition(const short_option_grammar &self) {
               
               using namespace boost::spirit::classic;
               // There's a boost::assertion too, so another explicit using
               // here:
               using boost::spirit::classic::assertion;
               
               // Assertions
               assertion<std::string> expect_root("expected object or array");
               assertion<std::string> expect_eoi("expected end of input");
               assertion<std::string> expect_objclose("expected ',' or ']'");
               
               // Grammar rules
               root 
                    =   boost::spirit::classic::epsilon_p[typename Context::a_object_s(self.c)] 
                    >>  expect_root(list_p(member, ch_p(',')))
                    >>  boost::spirit::classic::epsilon_p[typename Context::a_object_e(self.c)] 
                    >>  expect_eoi(end_p)
                    ;
               
               object 
                    =   ch_p('[')
                    >>  boost::spirit::classic::epsilon_p[typename Context::a_object_s(self.c)] 
                    >> (ch_p(']')
                        >>  boost::spirit::classic::epsilon_p[typename Context::a_object_e(self.c)] 
                        | (list_p(member, ch_p(','))
                           >> expect_objclose(ch_p(']'))
                           >>  boost::spirit::classic::epsilon_p[typename Context::a_object_e(self.c)]
                             )
                         )
                    ;

               member 
                    = symbol[typename Context::a_name(self.c)]
                    >> (ch_p(':')
                        >> value
                        | boost::spirit::classic::epsilon_p[typename Context::a_string_val(self.c)] )
                    ;
                
               value 
                    =   string[typename Context::a_string_val(self.c)] 
                    | (number | str_p("true") | "false" | "null")[typename Context::a_literal_val(self.c)]
                    | object
                    ;
                
               number 
                    =   !ch_p("-") >>
                    (ch_p("0") | (range_p(Ch('1'), Ch('9')) >> *digit_p)) >>
                    !(ch_p(".") >> +digit_p) >>
                    !(chset_p(boost::property_tree::detail::widen<Ch>("eE").c_str()) >>
                      !chset_p(boost::property_tree::detail::widen<Ch>("-+").c_str()) >>
                      +digit_p)
                    ;

               symbol 
                    =   lexeme_d[+character]
                    ;

               string
                    =   +(lexeme_d[confix_p('\"', *character, '\"')])
                    ;

               character
                    =   (anychar_p - "\\" - "\"" - ":" - "]" - ",")
                    [typename Context::a_char(self.c)]
                    ;


               // Debug
               BOOST_SPIRIT_DEBUG_RULE(root);
               BOOST_SPIRIT_DEBUG_RULE(object);
               BOOST_SPIRIT_DEBUG_RULE(member);
               BOOST_SPIRIT_DEBUG_RULE(value);
               BOOST_SPIRIT_DEBUG_RULE(symbol);
               BOOST_SPIRIT_DEBUG_RULE(string);
               BOOST_SPIRIT_DEBUG_RULE(number);
               BOOST_SPIRIT_DEBUG_RULE(character);

          }

          const boost::spirit::classic::rule<Scanner> &start() const {
               return root;
          }

     };

};

template<class It, class Ch>
unsigned long count_lines(It begin, It end)
{
     return static_cast<unsigned long>(std::count(begin, end, Ch('\n')) + 1);
}


template<class Ptree>
void read_short_option_format(std::basic_istream<typename Ptree::key_type::value_type> &stream, Ptree &pt) {
     
     using namespace boost::spirit::classic;
     typedef typename Ptree::key_type::value_type Ch;
     typedef typename std::vector<Ch>::iterator It;

     // Load data into vector
     std::vector<Ch> v(std::istreambuf_iterator<Ch>(stream.rdbuf()),
                       std::istreambuf_iterator<Ch>());
     if (!stream.good()) {
          boost::throw_exception(boost::property_tree::file_parser_error("read error", "", 0)); 
          std::exit(1);
     }
        
     // Prepare grammar
     short_option_grammar<Ptree> g;

     // Parse
     try {
          parse_info<It> pi = parse(v.begin(), v.end(), g, 
                                    space_p | comment_p("//") | comment_p("/*", "*/"));
          if (!pi.hit || !pi.full) {
               boost::throw_exception(parser_error<std::string, It>(v.begin(), "syntax error"));
               exit(1);
          }
     }
     catch (parser_error<std::string, It> &e)
     {
          boost::throw_exception(boost::property_tree::file_parser_error(e.descriptor, "", count_lines<It, Ch>(v.begin(), e.where)));
          exit(1);
     }

     // Swap grammar context root and pt
     pt.swap(g.c.root);

}


}}

#endif

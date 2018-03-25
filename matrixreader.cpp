#include "matrixreader.h"

MatrixReader::MatrixReader( std::string const & filename )
{
    std::vector<std::string> content;
    loadFormatFile( filename, content );
    parseFormatFile( content );
}

void MatrixReader::loadFormatFile( std::string const & filename, std::vector<std::string> & content )
{
    std::ifstream infile( filename );
    if ( !infile )
        throw std::invalid_argument("Can't open a matrix format file");

    std::string tmpString;

    while( getline(infile, tmpString) )
        content.push_back( tmpString );

    infile.close();
}

inline bool MatrixReader::isInteger(const std::string & s)
/*
 *strtol seems quite raw at first glance, so an explanation will make the code simpler to read :
  strtol will parse the string, stopping at the first character that cannot be considered part of an integer.
  If you provide p (as I did above), it sets p right at this first non-integer character.
  My reasoning is that if p is not set to the end of the string (the 0 character),
    then there is a non-integer character in the string s, meaning s is not a correct integer.
 The first tests are there to eliminate corner cases (leading spaces, empty string, etc.).
 */
{
   if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false ;

   char * p ;
   strtol(s.c_str(), &p, 10) ;

   return (*p == 0) ;
}

void MatrixReader::parseFormatFile( std::vector<std::string> & content )
{
    bool is_empty;
    std::string variable, value;

    for ( size_t k = 0; k < content.size(); k++ )
    {
        parseLine( content[k], is_empty, variable, value, k );
        if ( !is_empty )
        {
            if ( variable == "factor" ) factor = string_to_int(value, k );
            else if ( variable == "d_power" ) d_power = string_to_int(value, k);
            else if ( isInteger(variable) && isInteger(value) ) // если оба поля являются интами
                mtxFmt.emplace_back( std::stoi(variable), std::stoi(value) );
            else
                throw std::invalid_argument("Unknown variable: " + variable);

            value.clear();
            variable.clear();
        }
    }
}

void MatrixReader::parseLine( std::string & line, bool & is_empty, std::string & variable, std::string & value, int lineNumber )
{
    // изначально считаем, что строка пустая
    is_empty = true;

    // удаляем комментарии
    size_t pos = line.find('%');
    if ( pos != std::string::npos )
        line.erase( pos, line.length() - pos); // начиная откуда удалять, количество удаленных символов

    // если ничего кроме пробельных символов нет, то выходим
    std::string space_symbols = "\t \n";
    pos = line.find_first_not_of(space_symbols);
    if ( pos == std::string::npos )
        return;

    // ищем разделитель -- двоеточие
    pos = line.find(':');
    if ( pos != std::string::npos )
    {
        // если нашли, то строка не пустая
        is_empty = false;

        // смотрим на левую часть
        std::string lhs = line.substr(0, pos);

        // локализуем название перемнной
        size_t variable_start = lhs.find_first_not_of(space_symbols);
        size_t variable_end = lhs.find_last_not_of(space_symbols);

        variable = lhs.substr(variable_start, variable_end - variable_start + 1);
        //std::cout << "variable = " << variable << std::endl;

        // смотрим на правую часть
        std::string rhs = line.substr(pos + 1, line.length());

        size_t value_start = rhs.find_first_not_of(space_symbols);
        size_t value_end = rhs.find_last_not_of(space_symbols);
        value = rhs.substr(value_start, value_end - value_start + 1);
        //std::cout << "value = " << value << std::endl;
    }
    else
    {
        std::string stringLineNumber = std::to_string(lineNumber);
        throw std::invalid_argument("There is no assignment on line " + stringLineNumber);
    }
}

int MatrixReader::string_to_int( std::string const & value, size_t line )
{
    int result;
    if ( std::sscanf( value.c_str(), "%d", &result ) != 1 )
    {
        std::string lineNumberString = std::to_string( line );
        throw std::invalid_argument("Can't transform string to int in line " + lineNumberString );
    }

    return result;
}

void MatrixReader::showMatrixFormat()
{
    std::cout << "factor: " << factor << std::endl;
    std::cout << "d_power: " << d_power << std::endl;
    for ( size_t k = 0; k < mtxFmt.size(); k++ )
        std::cout << "diagonal: " << mtxFmt[k].first << "; value: " << mtxFmt[k].second << std::endl;
}

void MatrixReader::fillMatrix(Eigen::MatrixXd & m, double d)
{
    int size = m.rows();
    std::cout << "(MatrixReader) factor: " << factor << "; d_power: " << d_power << std::endl;

    for ( size_t k = 0; k < mtxFmt.size(); k++ )
    {
        int i = 0;
        int j = mtxFmt[k].first;
        double value = mtxFmt[k].second / factor * pow(d, d_power);

        for ( ; j < size; i++, j++ )
            m(i, j) = value;

        i = mtxFmt[k].first;
        j = 0;
        for ( ; i < size; i++, j++ )
            m(i, j) = value;
    }

}

void MatrixReader::resetFile(const std::string &filename)
{
    // чистим форматный вектор
    mtxFmt.clear();

    std::vector<std::string> content;
    loadFormatFile( filename, content );
    parseFormatFile( content );
}

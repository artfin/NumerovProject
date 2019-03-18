#include "matrixreader.h"

MatrixReader::MatrixReader( std::string const & filename ) : filename(filename)
{
    std::vector<std::string> content;
    loadFormatFile( );
    parseContent( );
}

void MatrixReader::loadFormatFile( )
{
    std::ifstream infile( filename );
    if ( !infile )
        throw std::invalid_argument("Can't open a matrix format file");

    std::string tmpString;

    while( getline(infile, tmpString) )
        content.push_back( tmpString );

    infile.close();
}

std::pair<int, int> MatrixReader::parseTwoIntegers(const std::string &s, size_t line)
{
   size_t pos = s.find(' ');
   std::string integer1 = s.substr(0, pos);
   std::string integer2 = s.substr(pos + 1, s.size());

   return std::make_pair( string_to_int(integer1, line),
                          string_to_int(integer2, line));
}


inline bool MatrixReader::isTwoIntegers(const std::string & s)
{
    size_t pos = s.find(' ');
    if ( pos == std::string::npos) {
        return false;
    }

    std::string integer1 = s.substr(0, pos);
    if ( !isInteger(integer1) ) {
        return false;
    }

    std::string integer2 = s.substr(pos+1, s.size());
    if ( !isInteger(integer2) ) {
        return false;
    }

    return true;
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

void MatrixReader::parseContent( )
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
            else if ( isTwoIntegers(variable) ) {
                //std::cout << "Two integers!" << std::endl;
                std::pair<int, int> integers = parseTwoIntegers(variable, k);

                extreme.insert( std::make_pair(integers, std::stoi(value)) );
            }
            else if ( isInteger(variable) && isInteger(value) ) // если оба поля являются интами
                mtxFmt.emplace_back( std::stoi(variable), std::stoi(value) );
            else
                throw std::invalid_argument("Unknown variable: " + variable);

            value.clear();
            variable.clear();
        }
    }
}

void MatrixReader::parseLine( std::string & line, bool & is_empty, std::string & variable, std::string & value, size_t lineNumber )
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
        //std::cout << "lhs: " << lhs << std::endl;

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

void MatrixReader::fillMatrix(Eigen::MatrixXd & m, double d)
{
    size_t size = m.rows();
    //std::cout << "(MatrixReader) factor: " << factor << "; d_power: " << d_power << std::endl;

    for ( auto const& p : mtxFmt )
    {
        size_t i = 0;
        size_t j = p.first;

        double value = p.second / factor * std::pow(d, d_power);

        for ( ; j < size; i++, j++ )
            m(i, j) = value;

        i = p.first;
        j = 0;
        for ( ; i < size; i++, j++ )
            m(i, j) = value;
    }

    for ( auto const & p : extreme ) {
        double value = p.second / factor * std::pow(d, d_power);

        int index1 = p.first.first;
        int index2 = p.first.second;

        m(index1, index2) = value;
        m(index2, index1) = value;

        m(size-1-index1, size-1-index2) = value;
        m(size-1-index2, size-1-index1) = value;
    }
}

void MatrixReader::resetFile(std::string const& filename)
{
    // чистим форматный вектор
    mtxFmt.clear();
    content.clear();
    extreme.clear();

    this->filename = filename;
    loadFormatFile();
    parseContent();
}

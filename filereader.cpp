#include "filereader.h"

FileReader::FileReader(std::string filename, Parameters * parameters)
    : filename(filename), parameters(parameters)
{
    std::cout << "FileReader constructor" << std::endl;

    std::vector<std::string> content;
    loadContent( content );
    parseContent( content );
}

void FileReader::loadContent( std::vector<std::string> & content )
{
    std::ifstream infile(filename);
    if (!infile)
    {
        throw std::invalid_argument("Can't open the file!");
    }

    std::string tmpString;

    while ( getline(infile, tmpString) )
        content.push_back( tmpString );

    infile.close();
}

void FileReader::parseContent( std::vector<std::string> & content )
{
    bool is_assignment; // есть ли присваивание переменной в текущей строке
    std::string variable, value;

    for ( size_t k = 0; k < content.size(); k++ )
    {
        parseLine( content[k], is_assignment, variable, value, k );

       // если присваивание происходит, то заполняем поля структуры
        if ( is_assignment )
        {
            if ( variable == "d") parameters->set_d( string_to_double(value, k) );
            else if ( variable == "mass" ) parameters->set_mass( string_to_double(value, k) );
            else if ( variable == "maxEnergy") parameters->set_maxEnergy( string_to_double(value, k) );
            else if ( variable == "epsilon") parameters->set_epsilon( string_to_double(value, k) );
            else if ( variable == "lowerBound" ) parameters->set_lowerBound( string_to_double(value, k) );
            else if ( variable == "upperBound" ) parameters->set_upperBound( string_to_double(value, k) );
            else
            {
                throw std::invalid_argument("Unknown variable: " + variable);
            }

            value.clear();
            variable.clear();
        }
    }
}

void FileReader::parseLine( std::string line, bool & is_assignment, std::string & variable, std::string & value, int lineNumber )
{
    is_assignment = false;

    // если нашли символ %, то удаляем все, что за ним -- это комментарий
    size_t pos = line.find('%');
    if ( pos != std::string::npos )
    {
        line.erase( pos, line.length() - pos );
    }

    std::string space_symbols = "\n \t";
    // ищем первый не пробельный символ
    // если дошли до конца, значит строка пустая, выходим из функции
    pos = line.find_first_not_of(space_symbols);
    if ( pos == std::string::npos )
        return;


    // если нашли знак равно, значит у нас здесь происходит присваивание переменной
    pos = line.find('=');
    if ( pos != std::string::npos )
    {
        is_assignment = true;

        // отрезаем левую часть перед '='
        std::string lhs = line.substr(0, pos);

        // ограждаем название переменной справа и слева
        size_t variable_start = lhs.find_first_not_of(space_symbols);
        size_t variable_end = lhs.find_last_not_of(space_symbols);
        variable = lhs.substr(variable_start, variable_end - variable_start + 1);
        //std::cout << "variable = " << variable << std::endl;

        // отрезаем правую часть
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

int FileReader::string_to_int( std::string const & value, size_t line )
{
    int result;
    if ( std::sscanf( value.c_str(), "%d", &result ) != 1 )
    {
        std::string lineNumberString = std::to_string( line );
        throw std::invalid_argument("Can't transform string to int in line " + lineNumberString );
    }

    return result;
}

double FileReader::string_to_double( std::string const & value, size_t line )
{
    double result;
    if ( std::sscanf( value.c_str(), "%lg", &result ) != 1 )
    {
        std::string lineNumberString = std::to_string( line );
        throw std::invalid_argument("Can't transform string to double in line " + lineNumberString );
    }

    return result;
}



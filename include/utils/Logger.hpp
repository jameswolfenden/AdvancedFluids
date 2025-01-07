#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>
#include <mutex>

namespace fluid
{

    enum class LogLevel
    {
        DEBUG,
        INFO,
        ERROR
    };

    class Logger
    {
    public:
        static void setLevel(LogLevel level);
        static void debug(const std::string &msg);
        static void info(const std::string &msg);
        static void error(const std::string &msg);

    private:
        static void log(LogLevel level, const std::string &msg);
        static std::mutex mutex_;
        static LogLevel current_level_;
    };

} // namespace fluid

#endif
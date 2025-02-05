#include "utils/Logger.hpp"
#include <iostream>

namespace fluid
{

    std::mutex Logger::mutex_;
    LogLevel Logger::current_level_ = LogLevel::INFO; // Uses the enum value of INFO

    void Logger::setLevel(LogLevel level)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        current_level_ = level;
    }

    void Logger::debug(const std::string &msg)
    {
        log(LogLevel::DEBUG, msg);
    }

    void Logger::info(const std::string &msg)
    {
        log(LogLevel::INFO, msg);
    }

    void Logger::error(const std::string &msg)
    {
        log(LogLevel::ERROR, msg);
    }

    void Logger::log(LogLevel level, const std::string &msg)
    {
        if (level < current_level_)
        {
            return; // Don't output messages e.g. if the level is set to ERROR and the message is INFO
        }

        std::lock_guard<std::mutex> lock(mutex_); // Lock the mutex to prevent multiple threads from writing to the console at the same time

        switch (level)
        {
        case LogLevel::DEBUG:
            std::cout << "DEBUG: " << msg << std::endl;
            break;
        case LogLevel::INFO:
            std::cout << msg << std::endl;
            break;
        case LogLevel::ERROR:
            std::cerr << "ERROR: " << msg << std::endl;
            break;
        }
    }
} // namespace fluid